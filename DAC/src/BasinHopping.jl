#include("../base/Atoms.jl")
#include("../base/FIRE.jl")
#include("../base/CNA.jl")

struct BasinHopper
	optimizer::Optimizer
	calculator::Calculator
	metC::MetC
	reseeder::Reseeder
	formula::Dict{String, Int64}
	boxLength::Float64
	vacuumAdd::Float64
	kT::Float64
	perturber::Function
	fmax::Float64
	rcut::Float64
	walltime::Float64
	recordingMode::String
	io::Tuple{IO, Channel}
	logIO::Tuple{IO, Channel}
	CNAIO::Tuple{IO, Channel}
	clusterVector::ClusterVector
	logResumeFile::Bool
	exitOnReseed::Bool
	calculateNormalCNA::Bool
	version::String
end

function logCNA(io::Tuple{IO, Channel}, ID::Int64, CNA::CNAProfile, energy::Float64)
	# block the file

	s = "$ID="
	for pair in CNA
		s *= "($(string(pair.first[1])),$(string(pair.first[2])),$(string(pair.first[3]))):$(string(pair.second));"
	end
	# unblock
	s *= "E$energy\n"

	print(io[1], s)

	return nothing
end


function logStep(io::Tuple{IO, Channel}, walkID::Int64, step::Int64, clusterID::Int64, energy::Float64, accepted::String)
	print(io[1], "WalkID $(string(walkID)), Step $(string(step)), ID $(string(clusterID)), energy $(string(energy)), accpeted $(string(accepted))\n")
end

function logStep(io::Tuple{IO, Channel}, step::Int64, clusterID::Int64, energy::Float64, accepted::String)
	print(io[1], "Step $(string(step)), ID $(string(clusterID)), energy $(string(energy)), accpeted $(string(accepted))\n")
end


function logStep(io::Tuple{IO, Channel}, step::Int64, clusterID::Int64, energy::Float64, accepted::String, sim::Float64)
	print(io[1], "Step " * string(step) * ", ID " * string(clusterID) * ", energy " * string(energy) * ", accepted " * string(accepted) * ", similarity, " * string(sim) * "\n")
end

function addToVector!(cluster::Union{Cluster, ClusterCompressed}, clusterVector::ClusterVector, dp::Int64)
	#= binary search =#
	energy = round(cluster.energy, digits=dp)
	#println(clusterVector.N[])
	presentClusterID = 0
	L = 1
	R = clusterVector.N[]+1
	while L < R

		m = trunc(Int, (R + L)/2)
		try
			if clusterVector.vec[m].energy == energy
				R = m
				startM = m
				while m < clusterVector.N[]&& clusterVector.vec[m].energy == energy
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
		catch err
			println("ADDTOVECTOR!   $(length(clusterVector.vec)) $m $R $L $((R+L)/2) $(Threads.threadid()) $err")
			throw(err)
			#vecprint = open("vecprint.txt", "w")
			#println(vecprint, "ADDTOVECTOR!   $(length(clusterVector.vec)) $m $R $L $((R+L)/2) $(Threads.threadid()) $err")
			#(vecprint, "ADDTOVECTOR2!  $(clusterVector.vec)")
			#close(vecprint)
			#exit()
			#throw(err)
		end
	end	
	
	#println("$energy, $presentClusterID, $R")
	# If cluster is new, add it to the vector.
	if presentClusterID == 0
		R = R == 0 ? 1 : R
		presentClusterID = Threads.atomic_add!(clusterVector.N, 1)
		presentClusterID += 1
		insert!(clusterVector.vec, R, ClusterCompressed(cluster.positions, round(cluster.energy, digits=dp), cluster.CNA, presentClusterID))
	end

	# return 0 if the cluster was added, otherwise return the index in the vector the cluster was found at.
	return presentClusterID
end

function optRun(_opt::PyObject, workhorse::Workhorse, fmax::Float64)
	opt = _opt(workhorse._py, logfile=nothing)
	opt.run(fmax=fmax)
end


function hop(bh::BasinHopper, steps::Int64, stepsAtomic::Threads.Atomic{Int64}, seed::Union{String, Cluster}, walkID::Int64, additionalInfo::Dict{String, Any}, start::DateTime, version::String)
	if version != "v1.2.1" || bh.version != "v1.2.1"
		println(bh.io[1], "The version number passed to the hop function or BasinHopper constructor does not match\nthe hard coded
			version number. Double check you are using the correct run script. This program will now terminate.")
		return 0
	end

	#start = now()

	# If needed, generate a random seed.
	if seed == "random"
		seed = generateRandomSeed(bh.formula, bh.boxLength, bh.vacuumAdd)
	end
	
	oldCluster = deepcopy(seed)


	setCalculator!(oldCluster, bh.calculator)
	optimize!(bh.optimizer, oldCluster, bh.fmax)
	if bh.calculateNormalCNA
		setCNAProfiles!(oldCluster, bh.rcut)
	else
		setCNAProfile!(oldCluster, bh.rcut)
	end



	# Specify variables from additionalInfo.
	step::Int64							= haskey(additionalInfo, "stepsCompleted")			? additionalInfo["stepsCompleted"] 			: 0
	setHopsToReseed!(bh.reseeder, 		  haskey(additionalInfo, "hopsToReseed")			? additionalInfo["hopsToReseed"] 			: getReseedPeriod(bh.reseeder))
	setReseedEnergyToBeat!(bh.reseeder,   haskey(additionalInfo, "reseedEnergyToBeat")		? additionalInfo["reseedEnergyToBeat"] 		: Inf)
	Emin::Float64						= haskey(additionalInfo, "Emin") 					? additionalInfo["Emin"] 					: Inf
	EminLocatedAt::Int64				= haskey(additionalInfo, "EminLocatedAt")			? additionalInfo["EminLocatedAt"] 			: 0
	targets								= haskey(additionalInfo, "targets") 				? additionalInfo["targets"] 				: Vector{Float64}()
	targetCNAs							= haskey(additionalInfo, "targetCNAs") 				? additionalInfo["targetCNAs"] 				: Vector{CNAProfile}()
	targetRounding::Int64				= haskey(additionalInfo, "targetRounding") 			? additionalInfo["targetRounding"] 			: 2
	targetsFound 						= haskey(additionalInfo, "targetsFound") 			? additionalInfo["targetsFound"] 			: falses(length(targets))
	targetsLocatedAt					= haskey(additionalInfo, "targetsLocatedAt") 		? additionalInfo["targetsLocatedAt"] 		: zeros(Int64, length(targets))

	# Boolean for checking CNAs of targets
	checkCNAsOfTarget = length(targetCNAs) > 0


	# If there is at least on etarget, exit on locating it/them.
	exitOnLocatingTargets = length(targets) > 0 

	newCluster = Cluster(bh.formula, getPositions(oldCluster), getCell(oldCluster))
	setCalculator!(newCluster, bh.calculator)

	# Used for determining when the Garbage collector should be run.
	# Reseeds cause the `step` iteration to desync with the modulus.
	iterations = 1


	while step < steps
		println("$walkID start of step $step"); flush(stdout)
		# Break loop if walltime exceeded.
		if (now() - start) / Hour(1) > bh.walltime
			print(bh.io[1], "\nwallTime exceeded. Ending Walk $(walkID).")
			flush(bh.io[1])
			break
		end

		Threads.atomic_add!(stepsAtomic, 1)
		
		step += 1
		stepLog = ""
		stepLog *= "\n================================\n"
		stepLog *= "Attempting step $step in Walk $(walkID) with Walker $(Threads.threadid())"


		iterations += 1

		newPos, pertrubString = bh.perturber(getPositions(oldCluster))
		setPositions!(newCluster, newPos)
		

		println("$walkID A of step $step"); flush(stdout)
		optimize!(bh.optimizer, newCluster, bh.fmax)
		while !isClusterCoherent(newCluster.positions, 2)
			newPos, pertrubString = bh.perturber(getPositions(oldCluster))
			setPositions!(newCluster, newPos)
			optimize!(bh.optimizer, newCluster, bh.fmax)
		end
		calculateEnergy!(newCluster, bh.calculator)
		if bh.calculateNormalCNA
			setCNAProfiles!(newCluster, bh.rcut)
		else
			setCNAProfile!(newCluster, bh.rcut)
		end
			
		println("$walkID B of step $step"); flush(stdout)

		# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
		# clusterID will be negative if is was already in the vector.
		Threads.lock(bh.clusterVector.lock) do 
			clusterID = addToVector!(newCluster, bh.clusterVector, 2)
		end

		println("$walkID C of step $step"); flush(stdout)

		if clusterID > 0
			logCNA(bh.CNAIO, clusterID, getCNA(newCluster), getEnergy(newCluster))
			stepLog *= "\nGenerated new cluster:\n\tID = $clusterID\n\tE = $(getEnergy(newCluster))"
		else
			stepLog *=  "\nRegenerated cluster:\n\tID = $clusterID\n\tE = $(getEnergy(newCluster))"
		end

		println("$walkID D of step $step"); flush(stdout)

		# Determine if hop to new structure is to be accepted.
		acceptHop, MetCString = getAcceptanceBoolean(bh.metC, oldCluster, newCluster)
		stepLog *= MetCString
		acceptStr = acceptHop ? "accepted." : "rejected."

		stepLog *= "\nThe current step has been " * acceptStr

		println("$walkID E of step $step"); flush(stdout)

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
			if bh.calculateNormalCNA
				setCNAProfiles!(oldCluster, bh.rcut)
			else
				setCNAProfile!(oldCluster, bh.rcut)
			end
				
		end

		println("$walkID F of step $step"); flush(stdout)

		# Log the move.
		logStep(bh.logIO, walkID, step, abs(clusterID), newCluster.energy, string(acceptHop))

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
						stepLog *=  "\nTarget $(targets[t]) has been located at step $step"
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
				stepLog *= "\nAll targets have been located.\nRun ending.\n"
				break
			end
		end

		println("$walkID G of step $step"); flush(stdout)

		# Check if time for reseed. Will not trigger if hopsToReseed is negative.
		if timeToReseed!(bh.reseeder)
			# if this walk should exit upon a reseed:
			if bh.exitOnReseed
				stepLog *=  "$(getReseedPeriod(bh.reseeder)) steps have occured since the last improvement. endning walk $(walkID).\n"
				print(bh.io[1], stepLog)
				flush(bh.io[1])
				return step
			end
			step += 1 #treat the reseed as an additional hop.
			stepLog *=  "$(getReseedPeriod(bh.reseeder)) steps have occured since the last improvement. reseeding.\n"
			
			# Generate a new seed (only update the positions of oldCluster).

			setPositions!(oldCluster, bh.reseeder.getReseedStructure(bh.reseeder.args...))
			optimize!(bh.optimizer, oldCluster, bh.fmax)
			calculateEnergy!(oldCluster, bh.calculator)
			if bh.calculateNormalCNA
				setCNAProfiles!(oldCluster, bh.rcut)
			else
				setCNAProfile!(oldCluster, bh.rcut)
			end

			# Check if the cluster is unique, add it to the vector of clusters, and update the CNA log.
			clusterID = addToVector!(oldCluster, bh.clusterVector, 2)
			if clusterID > 0
				logCNA(bh.CNAIO, clusterID, getCNA(oldCluster), getEnergy(oldCluster))
				stepLog *= "\nGenerated new cluster:\n\tID = $clusterID\n\tE = $(getEnergy(oldCluster))"
			else
				stepLog *= "\nRegenerated cluster:\n\tID = $clusterID\n\tE = $(getEnergy(oldCluster))"
			end

			# Log the step. 
			logStep(bh.logIO, walkID, step, abs(clusterID), oldCluster.energy, "reseed")
			
			# Reset hopsToReseed.
			hopsToReseed = getReseedPeriod(bh.reseeder)
		end

		println("$walkID H of step $step"); flush(stdout)

		# log the step
		print(bh.io[1], stepLog)
		flush(bh.io[1])

	end

	# Log resumption information.
	if bh.logResumeFile
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
	end

	GC.gc()
	
	return step
end
