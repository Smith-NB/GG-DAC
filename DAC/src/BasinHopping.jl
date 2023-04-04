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

function hop(bh::BasinHopper, steps::Int64, seed::Union{String, Cluster}, additionalInfo::Dict{String, Union{Number, Vector{Float64}}})

	# If needed, generate a random seed.
	if seed == "random"
		seed = generateRandomSeed(bh.formula, bh.boxLength, bh.vacuumAdd)
	end
	
	# set some parameters and calculate energy and CNA profile of seed

	oldCluster = seed

	setPositions!(bh.workhorse, getPositions(seed))
	setCell!(bh.workhorse, getCell(seed))
	bh.workhorseOpt.run(fmax=bh.fmax)
	#optimize!(bh.optimizer, bh.workhorse, bh.fmax)

	setPositions!(oldCluster, getPositions(bh.workhorse))
	setCNAProfile!(oldCluster, bh.rcut)


	# specify remaining hops to reseed (the value stored in additionalInfo, if present, otherwise value from bh).
	hopsToReseed = haskey(additionalInfo, "hopsToReseed") ? additionalInfo["hopsToReseed"] : bh.reseedPeriod
	Emin = haskey(additionalInfo, "Emin") ? additionalInfo["Emin"] : Inf
	targets = haskey(additionalInfo, "targets") ? additionalInfo["targets"] : Vector{Float64}()
	targetRounding = haskey(additionalInfo, "targetRounding") ? additionalInfo["targetRounding"] : 2
	exitOnLocatingTargets = length(targets) > 0
	targetFound = falses(length(targets))
	hopsTargetsLocatedAt = zeros(length(targets))

	newCluster = Cluster(bh.formula, getPositions(oldCluster), getCell(oldCluster))

	# perform n steps
	for step in 1:steps
		println(bh.io[1], "\n================================\n")
		println(bh.io[1], "Attempting step ", step)

		#perturb cluster and recalculate energies and CNA profile
		#println("time to perturbCluster")
		#@time setPositions!(newCluster, perturbCluster(oldCluster.positions, bh.dr))
		#println("time to optimize!")
		#@time n = optimize!(bh.optimizer, newCluster, bh.fmax)
		#println(n, " calc calls")
		#calculateEnergy!(newCluster, newCluster.calculator)
		#println("time to setCNAProfile!")
		#@time setCNAProfile!(newCluster, bh.rcut)
		#println("time to addToVector!")
		#@time clusterIsUnique = addToVector!(newCluster, bh.clusterVector, 2)
		if step % 30 == 0
			GC.gc()
		end

		setPositions!(bh.workhorse, perturbCluster(getPositions(oldCluster), bh.dr))			
		ncalls = bh.workhorseOpt.run(fmax=bh.fmax)

		setPositions!(newCluster, getPositions(bh.workhorse))
		setEnergy!(newCluster, getEnergy!(bh.workhorse))
		setCNAProfile!(newCluster, bh.rcut)
		clusterIsUnique = addToVector!(newCluster, bh.clusterVector, 2)

		if clusterIsUnique
			logCNA(bh.CNAIO, step, getCNA(newCluster))
		end

		print(bh.io[1], "\nGenerated new cluster, E = ", getEnergy(newCluster))
		print(bh.io[1], "\n$ncalls made to py")

		# determine if hop to new structure is to be accepted
		acceptHop = getAcceptanceBoolean(bh.metC, oldCluster, newCluster)
		acceptStr = acceptHop ? "accepted." : "rejected."
		print(bh.io[1], "\nThe current step has been " * acceptStr)



		# if accepted, update
		if acceptHop
			# if new LES found, update Emin and reset hopsToReseed
			if getEnergy(newCluster) < Emin
				Emin = getEnergy(newCluster)
				hopsToReseed = bh.reseedPeriod
			end

			# update oldCluster to newCluster
			setPositions!(oldCluster, getPositions(newCluster))
			setEnergy!(oldCluster, getEnergy(newCluster))
			setCNAProfile!(oldCluster, getCNA(newCluster))
			
		else # if hop was rejected
			hopsToReseed -= 1
		end

		# Check if the new cluster is a target. Exit if all targets found.
		if acceptHop && exitOnLocatingTargets
			allTargetsFound = true
			for t in 1:length(targets)
				if round(getEnergy(oldCluster), digits=targetRounding) == targets[t]
					print(bh.io[1], "\nTarget ", targets[t], " has been located at step ", step)
					targetFound[t] = true
					hopsTargetsLocatedAt[t] = step
				end
				if !targetFound[t]
					allTargetsFound = false
				end
			end

			# all targets found, exit.
			if allTargetsFound
				print(bh.io[1], "\nAll targets have been located.\nRun ending.\n")
				return nothing
			end
		end

		logStep(bh.logIO, step, newCluster.energy, acceptHop)

		# Check if time for reseed. Will not trigger if hopsToReseed is negative.
		if hopsToReseed == 0
			print(bh.io[1], bh.reseedPeriod, " steps have occured since the last improvement. reseeding.")
			
			# generate a new seed (only update the positions of oldCluster)
			oldCluster = generateRandomSeed(bh.formula, bh.boxLength, bh.vacuumAdd, true)

			setPositions!(bh.workhorse, getPositions(seed))
			bh.workhorseOpt.run(fmax=bh.fmax)
			setPositions!(oldCluster, getPositions(bh.workhorse))
			setEnergy!(oldCluster, getEnergy!(bh.workhorse))
			setCNAProfile!(oldCluster, bh.rcut)

			clusterIsUnique = addToVector!(oldCluster, bh.clusterVector, 2)

			if clusterIsUnique
				logCNA(bh.CNAIO, step, getCNA(oldCluster))
			end

			# reset hopsToReseed
			hopsToReseed = bh.reseedPeriod
		end

	end

	
end
