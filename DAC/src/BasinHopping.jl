#include("../base/Atoms.jl")
#include("../base/FIRE.jl")
#include("../base/CNA.jl")

struct BasinHopper
	optimizer::Optimizer
	calculator::Calculator
	metC::MetC
	formula::Dict{String, Int64}
	boxLength::Number
	vacuumAdd::Number
	kT::Number
	dr::Number
	reseedPeriod::Int64
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
end

function logCNA(io::Tuple{IO, Channel}, ID::Int64, CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}})
	# create the string to log from the given CNA and cluster ID
	s = string(ID) * "="
	for pair in CNA
		s *= string(pair.first[1]) * "," * string(pair.first[2]) * "," * string(pair.first[3]) * ":" * string(pair.second) * ";"
	end
	# create break line
	s *= "\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end

	return nothing
end


function logStep(io::Tuple{IO, Channel}, step::Int64, energy::Float64, accepted::Bool)
	s = "Step " * string(step) * ", energy " * string(energy) * ", accepted " * string(accepted) * "\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
end

function logStep(io::Tuple{IO, Channel}, step::Int64, energy::Float64, accepted::Bool, sim::Float64)
	s = file, "Step " * string(step) * ", energy " * string(energy) * ", accepted " * string(accepted) * ", similarity, " * string(sim) * "\n"

	# add the string to the file channel. this line will block if the channel is full, i.e. by another thread
	put!(io[2], s)

	# take the string in the channel and print it to the output file.
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
end

function addToVector!(cluster::Cluster, clusterVector::ClusterVector, dp::Int64)
	#= binary search =#
	energy = round(cluster.energy, digits=dp)
	clusterAlreadyPresent = false
	L = 1
	R = clusterVector.N+1

	while L < R
		m = trunc(Int, (R + L)/2)

		if clusterVector.vec[m].energy == energy
			R = m
			startM = m
			while m < clusterVector.N && clusterVector.vec[m].energy == energy
				if clusterVector.vec[m].CNA == cluster.CNA
					clusterAlreadyPresent = true
					break
				end
				m += 1
			end

			m = startM - 1
			while m > 0 && clusterVector.vec[m].energy == energy
				if clusterVector.vec[m].CNA == cluster.CNA
					clusterAlreadyPresent = true
					break
				end
				m -= 1
			end
			#
			break

		elseif clusterVector.vec[m].energy < energy
			L = m + 1
		else
			R = m
		end
	end	

	if !clusterAlreadyPresent
		R = R == 0 ? 1 : R
		insert!(clusterVector.vec, R, ClusterCompressed(cluster.positions, round(cluster.energy, digits=dp), cluster.CNA))
		clusterVector.N += 1
	end

	# return true if an insertion to the vector was made
	return !clusterAlreadyPresent

end

function optRun(_opt::PyObject, workhorse::Workhorse, fmax::Float64)
	opt = _opt(workhorse._py, logfile=nothing)
	opt.run(fmax=fmax)
end

function logResumptionInfo()


end

function hop(bh::BasinHopper, steps::Int64, seed::Union{String, Cluster}, additionalInfo::Dict{String, Union{Number, Vector{Float64}}})

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
	step  					= haskey(additionalInfo, "stepsCompleted")			? additionalInfo["stepsCompleted"] 			: 0
	hopsToReseed  			= haskey(additionalInfo, "hopsToReseed")			? additionalInfo["hopsToReseed"] 			: bh.reseedPeriod
	reseedEnergyToBeat		= haskey(additionalInfo, "reseedEnergyToBeat")		? additionalInfo["reseedEnergyToBeat"] 		: Inf
	Emin  					= haskey(additionalInfo, "Emin") 					? additionalInfo["Emin"] 					: Inf
	EminLocatedAt			= haskey(additionalInfo, "EminLocatedAt")			? additionalInfo["EminLocatedAt"] 			: 0
	targets					= haskey(additionalInfo, "targets") 				? additionalInfo["targets"] 				: Vector{Float64}()
	targetRounding 	 		= haskey(additionalInfo, "targetRounding") 			? additionalInfo["targetRounding"] 			: 2
	targetsFound 			= haskey(additionalInfo, "targetsFound") 			? additionalInfo["targetsFound"] 			: falses(length(targets))
	targetsLocatedAt		= haskey(additionalInfo, "targetsLocatedAt") 		? additionalInfo["targetsLocatedAt"] 		: zeros(length(targets))

	# If there is at least on etarget, exit on locating it/them.
	exitOnLocatingTargets = length(targets) > 0 

	newCluster = Cluster(bh.formula, getPositions(oldCluster), getCell(oldCluster))

	
	while step < steps
		step += 1
		println(bh.io[1], "\n================================\n")
		println(bh.io[1], "Attempting step ", step)

		#= Manual Garbage collection (gross). When using asap3 as the workhorse, something goes wrong such that the refcount
		    of an object exceeds 100. This causes an assertion in C++ to fail. =#
		if step % 30 == 0
			GC.gc()
		end

		# Perturb and optimize with the workhorse.
		setPositions!(bh.workhorse, perturbCluster(getPositions(oldCluster), bh.dr))			
		bh.workhorseOpt.run(fmax=bh.fmax)

		# Update newCluster with optimized workhorse.
		setPositions!(newCluster, getPositions(bh.workhorse))
		setEnergy!(newCluster, getEnergy!(bh.workhorse))
		setCNAProfile!(newCluster, bh.rcut)

		# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
		clusterIsUnique = addToVector!(newCluster, bh.clusterVector, 2)
		if clusterIsUnique
			logCNA(bh.CNAIO, step, getCNA(newCluster))
		end

		print(bh.io[1], "\nGenerated new cluster, E = ", getEnergy(newCluster))

		# Determine if hop to new structure is to be accepted.
		acceptHop = getAcceptanceBoolean(bh.metC, oldCluster, newCluster)
		acceptStr = acceptHop ? "accepted." : "rejected."

		print(bh.io[1], "\nThe current step has been " * acceptStr)

		# If accepted, update.
		if acceptHop
			# If new LES found, update Emin and reset hopsToReseed.
			if getEnergy(newCluster) < Emin
				Emin = getEnergy(newCluster)
			end

			if getEnergy(newCluster) < reseedEnergyToBeat
				hopsToReseed = bh.reseedPeriod
			end

			# Update oldCluster to newCluster.
			setPositions!(oldCluster, getPositions(newCluster))
			setEnergy!(oldCluster, getEnergy(newCluster))
			setCNAProfile!(oldCluster, getCNA(newCluster))
			
		else # If hop was rejected.
			hopsToReseed -= 1
		end

		# Log the move.
		logStep(bh.logIO, step, newCluster.energy, acceptHop)

		# Check if the new cluster is a target. Exit if all targets found.
		if acceptHop && exitOnLocatingTargets
			allTargetsFound = true
			for t in 1:length(targets)
				if round(getEnergy(oldCluster), digits=targetRounding) == targets[t]
					print(bh.io[1], "\nTarget ", targets[t], " has been located at step ", step)
					targetsFound[t] = true
					targetsLocatedAt[t] = step
				end
				if !targetsFound[t]
					allTargetsFound = false
				end
			end

			# All targets found, exit.
			if allTargetsFound
				print(bh.io[1], "\nAll targets have been located.\nRun ending.\n")
				return nothing
			end
		end

		# Check if time for reseed. Will not trigger if hopsToReseed is negative.
		if hopsToReseed == 0

			print(bh.io[1], bh.reseedPeriod, " steps have occured since the last improvement. reseeding.\n")
			
			# Generate a new seed (only update the positions of oldCluster).
			setPositions!(bh.workhorse, generateRandomSeed(bh.formula, bh.boxLength, bh.vacuumAdd, true))
			bh.workhorseOpt.run(fmax=bh.fmax)
			setPositions!(oldCluster, getPositions(bh.workhorse))
			setEnergy!(oldCluster, getEnergy!(bh.workhorse))
			setCNAProfile!(oldCluster, bh.rcut)

			# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
			clusterIsUnique = addToVector!(oldCluster, bh.clusterVector, 2)
			if clusterIsUnique
				logCNA(bh.CNAIO, step, getCNA(oldCluster))
			end

			# Reset hopsToReseed.
			hopsToReseed = bh.reseedPeriod
			step += 1 #treat the reseed as an additional hop.
		end

	end

	# Log resumption information.
	resumeFile = open("informationForResuming.txt", "w")
	print(resumeFile, "stepsCompleted:$step\n")
	print(resumeFile, "Emin:$Emin\n")
	print(resumeFile, "EminLocatedAt:$EminLocatedAt\n")
	print(resumeFile, "reseedEnergyToBeat:$reseedEnergyToBeat\n")
	print(resumeFile, "hopsToReseed:$hopsToReseed\n")
	print(resumeFile, "targets:$targets\n")
	print(resumeFile, "targetsFound:$targetsFound\n")
	print(resumeFile, "targetsLocatedAt:$targetsLocatedAt\n")
	close(resumeFile)
	
end
