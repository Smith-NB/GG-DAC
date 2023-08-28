#include("../base/Atoms.jl")
#include("../base/FIRE.jl")
#include("../base/CNA.jl")

struct BasinHopper
	optimizer::Optimizer
	calculator::Calculator
	metC::MetC
	reseeder::Reseeder
	formula::Dict{String, Int64}
	boxLength::Number
	vacuumAdd::Number
	kT::Number
	perturber::Function
	fmax::Number
	rcut::Number
	walltime::Number
	recordingMode::String
	io::Tuple{IO, Channel}
	logIO::Tuple{IO, Channel}
	CNAIO::Tuple{IO, Channel}
	clusterVector::ClusterVector
	workhorse::Workhorse
	workhorseOpt::PyObject
	version::String
end

function logCNA(io::Tuple{IO, Channel}, ID::Int64, CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}, energy::Float64)
	# create the string to log from the given CNA and cluster ID
	s = string(ID) * "="
	for pair in CNA
		s *= string(pair.first[1]) * "," * string(pair.first[2]) * "," * string(pair.first[3]) * ":" * string(pair.second) * ";"
	end
	# create break line
	s *= "E$energy\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end

	return nothing
end


function logStep(io::Tuple{IO, Channel}, step::Int64, clusterID::Int64, energy::Float64, accepted::String)
	s = "Step " * string(step) * ", ID " * string(clusterID) * ", energy " * string(energy) * ", accepted " * string(accepted)  * "\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
end

function logStep(io::Tuple{IO, Channel}, step::Int64, clusterID::Int64, energy::Float64, accepted::String, sim::Float64)
	s = "Step " * string(step) * ", ID " * string(clusterID) * ", energy " * string(energy) * ", accepted " * string(accepted) * ", similarity, " * string(sim) * "\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
end

function addToVector!(cluster::Union{Cluster, ClusterCompressed}, clusterVector::ClusterVector, dp::Int64)
	#= binary search =#
	energy = round(cluster.energy, digits=dp)
	presentClusterID = 0
	L = 1
	R = clusterVector.N+1
	while L < R
		m = trunc(Int, (R + L)/2)

		if clusterVector.vec[m].energy == energy
			R = m
			startM = m
			while m < clusterVector.N && clusterVector.vec[m].energy == energy
				if clusterVector.vec[m].CNA == cluster.CNA
					presentClusterID = -clusterVector.vec[m].ID # ClusterId will be negative if already present in Vector
					break
				end
				m += 1
			end

			m = startM - 1
			while m > 0 && clusterVector.vec[m].energy == energy
				if clusterVector.vec[m].CNA == cluster.CNA
					presentClusterID = -clusterVector.vec[m].ID # ClusterId will be negative if already present in Vector
					break
				end
				m -= 1
			end

			break

		elseif clusterVector.vec[m].energy < energy
			L = m + 1
		else
			R = m
		end
	end	
	#println("$energy, $presentClusterID, $R")
	# If cluster is new, add it to the vector.
	if presentClusterID == 0
		R = R == 0 ? 1 : R
		clusterVector.N += 1
		presentClusterID = clusterVector.N # Set cluster ID to be the new ID of the cluster
		insert!(clusterVector.vec, R, ClusterCompressed(cluster.positions, round(cluster.energy, digits=dp), cluster.CNA, presentClusterID))
	end

	# return 0 if the cluster was added, otherwise return the index in the vector the cluster was found at.
	return presentClusterID
end

function optRun(_opt::PyObject, workhorse::Workhorse, fmax::Float64)
	opt = _opt(workhorse._py, logfile=nothing)
	opt.run(fmax=fmax)
end


function hop(bh::BasinHopper, steps::Int64, seed::Union{String, Cluster}, additionalInfo::Dict{String, Any}, version::String)

	if version != "v1.2.1" || bh.version != "v1.2.1"
		println(bh.io[1], "The version number passed to the hop function or BasinHopper constructor does not match\nthe hard coded
			version number. Double check you are using the correct run script. This program will now terminate.")
		return nothing
	end

	start = now()

	# If needed, generate a random seed.
	if seed == "random"
		seed = generateRandomSeed(bh.formula, bh.boxLength, bh.vacuumAdd)
	end
	
	oldCluster = seed

	# Set the positions of the workhorse.
	setPositions!(bh.workhorse, getPositions(seed))
	setCell!(bh.workhorse, getCell(seed))
	bh.workhorseOpt.run(fmax=bh.fmax)

	# Update oldCluster with the optimized positions.
	setPositions!(oldCluster, getPositions(bh.workhorse))
	setEnergy!(oldCluster, getEnergy!(bh.workhorse))
	setCNAProfile!(oldCluster, bh.rcut)


	# Specify variables from additionalInfo.
	step  								= haskey(additionalInfo, "stepsCompleted")			? additionalInfo["stepsCompleted"] 			: 0
	setHopsToReseed!(bh.reseeder, 		  haskey(additionalInfo, "hopsToReseed")			? additionalInfo["hopsToReseed"] 			: getReseedPeriod(bh.reseeder))
	setReseedEnergyToBeat!(bh.reseeder,   haskey(additionalInfo, "reseedEnergyToBeat")		? additionalInfo["reseedEnergyToBeat"] 		: Inf)
	Emin  								= haskey(additionalInfo, "Emin") 					? additionalInfo["Emin"] 					: Inf
	EminLocatedAt						= haskey(additionalInfo, "EminLocatedAt")			? additionalInfo["EminLocatedAt"] 			: 0
	targets								= haskey(additionalInfo, "targets") 				? additionalInfo["targets"] 				: Vector{Float64}()
	targetCNAs							= haskey(additionalInfo, "targetCNAs") 				? additionalInfo["targetCNAs"] 				: Vector{CNAProfile}()
	targetRounding 	 					= haskey(additionalInfo, "targetRounding") 			? additionalInfo["targetRounding"] 			: 2
	targetsFound 						= haskey(additionalInfo, "targetsFound") 			? additionalInfo["targetsFound"] 			: falses(length(targets))
	targetsLocatedAt					= haskey(additionalInfo, "targetsLocatedAt") 		? additionalInfo["targetsLocatedAt"] 		: zeros(Int64, length(targets))

	# Boolean for checking CNAs of targets
	checkCNAsOfTarget = length(targetCNAs) > 0


	# If there is at least on etarget, exit on locating it/them.
	exitOnLocatingTargets = length(targets) > 0 

	newCluster = Cluster(bh.formula, getPositions(oldCluster), getCell(oldCluster))

	clusterVectorChannel = Channel{Cluster}(0) # unbuffered channel. put! blocks until matching !take is called.

	# Used for determining when the Garbage collector should be run.
	# Reseeds cause the `step` iteration to desync with the modulus.
	iterations = 1

	clusterVectorLock = ReentrantLock()

	while step < steps
		
		# Break loop if walltime exceeded.
		if (now() - start) / Hour(1) > bh.walltime
			break
		end

		step += 1
		stepLog = ""
		stepLog += "\n================================\n")
		stepLog += "Attempting step $step with Walker $(Threads.nthreads())"

		#= Manual Garbage collection (gross). When using asap3 as the workhorse, something goes wrong such that the refcount
		    of an object exceeds 100. This causes an assertion in C++ to fail. =#
		if iterations % 30 == 0
			GC.gc()
		end

		iterations += 1

		# Perturb and optimize with the workhorse.
		setPositions!(bh.workhorse, bh.perturber(getPositions(oldCluster)))
		stepLog += "\nOldEnergy = $(getEnergy(oldCluster))"
		bh.workhorseOpt.run(fmax=bh.fmax)

		# Update newCluster with optimized workhorse.
		setPositions!(newCluster, getPositions(bh.workhorse))
		setEnergy!(newCluster, getEnergy!(bh.workhorse))
		setCNAProfile!(newCluster, bh.rcut)

		# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
		# clusterID will be negative if is was already in the vector.
		begin
			lock(clusterVectorLock)
			try
				put!(clusterVectorChannel, newCluster)
				clusterID = addToVector!(take!(clusterVectorChannel, bh.clusterVector, 2))
			finally
				unlock(clusterVectorLock)
			end
		end

		if clusterID > 0
			logCNA(bh.CNAIO, clusterID, getCNA(newCluster), getEnergy(newCluster))
			stepLog += "\nGenerated new cluster:\n\tID = $clusterID\n\tE = $(getEnergy(newCluster))"
		else
			stepLog +=  "\nRegenerated cluster:\n\tID = $clusterID\n\tE = $(getEnergy(newCluster))"
		end

		# Determine if hop to new structure is to be accepted.
		acceptHop, MetCString = getAcceptanceBoolean(bh.metC, oldCluster, newCluster)
		stepLog += MetCString
		acceptStr = acceptHop ? "accepted." : "rejected."

		stepLog += "\nThe current step has been " * acceptStr

		# Decrease number of hops until next reseed.
		updateHopsToReseed!(bh.reseeder)

		# If accepted, update.
		if acceptHop
			# If new LES found, update Emin and reset hopsToReseed.
			if getEnergy(newCluster) < Emin
				Emin = getEnergy(newCluster)
			end

			# Check if a new LES has been found since the last reseed.
			checkNewlyAcceptedStructure!(bh.reseeder, newCluster)

			# Update oldCluster to newCluster.
			setPositions!(oldCluster, getPositions(newCluster))
			setEnergy!(oldCluster, getEnergy(newCluster))
			setCNAProfile!(oldCluster, getCNA(newCluster))
			
		end

		# Log the move.
		logStep(bh.logIO, step, abs(clusterID), newCluster.energy, string(acceptHop))

		# Check if the new cluster is a target. Exit if all targets found.
		if acceptHop && exitOnLocatingTargets
			allTargetsFound = true
			for t in 1:length(targets)
				if round(getEnergy(oldCluster), digits=targetRounding) == targets[t]
					
					found = false
					# If CNA checking enabled, ensure CNAs match
					if checkCNAsOfTarget && getCNASimilarity(targetCNAs[t], getCNAProfile(oldCluster)) == 1.0
						found = true
					elseif !checkCNAsOfTarget # If CNA checking disabled, assume target is found.
						found = true
					end

					if found
						stepLog +=  "\nTarget $(targets[t]) has been located at step $step"
						targetsFound[t] = true
						targetsLocatedAt[t] = step
					end
				end

				if !targetsFound[t]
					allTargetsFound = false
				end
			end

			# All targets found, exit.
			if allTargetsFound
				stepLog += "\nAll targets have been located.\nRun ending.\n")
				break
			end
		end

		# Check if time for reseed. Will not trigger if hopsToReseed is negative.
		if timeToReseed!(bh.reseeder)
			step += 1 #treat the reseed as an additional hop.
			stepLog +=  "$(getReseedPeriod(bh.reseeder)) steps have occured since the last improvement. reseeding.\n")
			
			# Generate a new seed (only update the positions of oldCluster).
			setPositions!(bh.workhorse, bh.reseeder.getReseedStructure(bh.reseeder.args...))
			bh.workhorseOpt.run(fmax=bh.fmax)
			setPositions!(oldCluster, getPositions(bh.workhorse))
			setEnergy!(oldCluster, getEnergy!(bh.workhorse))
			setCNAProfile!(oldCluster, bh.rcut)

			# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
			clusterID = addToVector!(oldCluster, bh.clusterVector, 2)
			if clusterID > 0
				logCNA(bh.CNAIO, clusterID, getCNA(oldCluster), getEnergy(oldCluster))
				stepLog += "\nGenerated new cluster:\n\tID = $clusterID\n\tE = $(getEnergy(oldCluster))")
			else
				stepLog += "\nRegenerated cluster:\n\tID = $clusterID\n\tE = $(getEnergy(oldCluster))")
			end

			# Log the step. 
			logStep(bh.logIO, step, abs(clusterID), oldCluster.energy, "reseed")
			
			# Reset hopsToReseed.
			hopsToReseed = getReseedPeriod(bh.reseeder)
		end

		# log the step
		put!(bh.io[2], stepLog)
		print(bh.io[1], take!(bh.io[2]))

	end

	# Log resumption information.
	resumeFile = open("informationForResuming.txt", "w")
	print(resumeFile, "stepsCompleted:$step\n")
	print(resumeFile, "Emin:$Emin\n")
	print(resumeFile, "EminLocatedAt:$EminLocatedAt\n")
	print(resumeFile, "reseedEnergyToBeat:$(getReseedEnergyToBeat(bh.reseeder))\n")
	print(resumeFile, "hopsToReseed:$(getHopsToReseed(bh.reseeder))\n")
	print(resumeFile, "targets:$targets\n")
	print(resumeFile, "targetsFound:$targetsFound\n")
	print(resumeFile, "targetsLocatedAt:$targetsLocatedAt\n")
	close(resumeFile)
	
	return nothing
end
