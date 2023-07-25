abstract type MetC end


#=============================================================================#
#=================================EnergyMetC==================================#
#=============================================================================#

struct EnergyMetC <: MetC
	kT::Float64
	io::Tuple{IO, Channel}
end


"""
	getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on the EnergyMetC.
"""
function getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)
	
	if newCluster.energy < oldCluster.energy
		return true
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	
	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(MetC.io[2], "\nChance to accept = " * string(probability))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end
	
	
	accept = probability > rand()
	return accept

end


#=============================================================================#
#============================getAcceptanceBoolean=============================#
#=============================================================================#

struct EnergyWithReferenceRestrictionMetC <: MetC
	kT::Float64
	refCNA::CNAProfile
	sigmaCut::Float64
	io::Tuple{IO, Channel}
end

"""
	getAcceptanceBoolean(MetC::EnergyWithReferenceRestrictionMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on the EnergyMetC, with the added restriction that if the newCluster is less
	similar to some reference than `sigmaCut`, the hop is rejected.
"""
function getAcceptanceBoolean(MetC::EnergyWithReferenceRestrictionMetC, oldCluster::Cluster, newCluster::Cluster)
	s = getCNASimilarity(MetC.refCNA, getCNAProfile(newCluster))
	if s < MetC.sigmaCut
		put!(MetC.io[2], "\nChance to accept = 0\nSimilarity to ref: $s")
		while isready(MetC.io[2])
			print(MetC.io[1], take!(MetC.io[2]))
		end
		return false
	end

	if newCluster.energy < oldCluster.energy
		return true
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	
	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(MetC.io[2], "\nChance to accept = " * string(probability))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end
	
	accept = probability > rand()
	return accept

end


#=============================================================================#
#===========================EnergyAndStructureMetC============================#
#=============================================================================#

struct EnergyAndStructureMetC <: MetC
	cE::Float64
	cSCM::Float64
	eFunction::Function
	simFunction::Function
	autoAcceptDownhillMove::Bool
	io::Tuple{IO, Channel}
end

"""
	getAcceptanceBoolean(MetC::EnergyAndStructureMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on the EnergyAndStructureMetC.
"""
function getAcceptanceBoolean(MetC::EnergyAndStructureMetC, oldCluster::Cluster, newCluster::Cluster)

	if MetC.autoAcceptDownhillMove && getEnergy(newCluster) < getEnergy(oldCluster)
		return true
	end

	eScore = MetC.eFunction(getEnergy(oldCluster), getEnergy(newCluster)) * MetC.cE
	simScore = MetC.simFunction(getCNASimilarity(getCNA(oldCluster), getCNA(newCluster))) * MetC.cSCM

	put!(MetC.io[2], "\nChance to accept = " * string(eScore + simScore))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end

	accept = (eScore + simScore) > rand()
	return accept

end	


#=============================================================================#
#=================================HISTOMetC===================================#
#=============================================================================#

# Accounts for accepted hops only

mutable struct HISTOMetC <: MetC
	kT::Float64
	w::Float64
	delta::Float64
	resetPeriod::Float64
	refCNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOMetC, oldCluster::Cluster, newCluster::Cluster)

	if newCluster.energy < oldCluster.energy
		return true
	end

	updateHist = false
	binNew = nothing
	put!(MetC.io[2], "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)")
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			put!(MetC.io[2], "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)")
		end
		MetC.timeElapsed += 1

		hScore = 0
	else
		histSum = sum(MetC.hist) # used to normalise bar heights
		if histSum == 0 # avoid division by zero
			histSum = 1
		end

		simOld = getCNASimilarity(getCNA(oldCluster), MetC.refCNA)
		binOld = simOld == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simOld/MetC.delta) + 1 # get bin of old Cluster
		hOld   = MetC.hist[binOld]/histSum # get height (normalised) of bars
		
		put!(MetC.io[2], "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)")

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		put!(MetC.io[2], "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)")

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = exp((oldCluster.energy - newCluster.energy + hScore) / MetC.kT)	

	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(MetC.io[2], "\nhScore = $(hScore)\nChance to accept = " * string(probability))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end
	
	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	end

	return accept

end


#=============================================================================#
#================================HISTOAttMetC=================================#
#=============================================================================#

# Accounts for all attempted hops

mutable struct HISTOAttMetC <: MetC
	kT::Float64
	w::Float64
	delta::Float64
	resetPeriod::Float64
	refCNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOAttMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOAttMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOAttMetC, oldCluster::Cluster, newCluster::Cluster)

	if newCluster.energy < oldCluster.energy
		return true
	end

	updateHist = false
	binNew = nothing
	put!(MetC.io[2], "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)")
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			put!(MetC.io[2], "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)")
		end
		MetC.timeElapsed += 1

		hScore = 0
	else
		histSum = sum(MetC.hist) # used to normalise bar heights
		if histSum == 0 # avoid division by zero
			histSum = 1
		end

		simOld = getCNASimilarity(getCNA(oldCluster), MetC.refCNA)
		binOld = simOld == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simOld/MetC.delta) + 1 # get bin of old Cluster
		hOld   = MetC.hist[binOld]/histSum # get height (normalised) of bars
		
		put!(MetC.io[2], "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)")

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		put!(MetC.io[2], "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)")

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = exp((oldCluster.energy - newCluster.energy + hScore) / MetC.kT)	

	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(MetC.io[2], "\nhScore = $(hScore)\nChance to accept = " * string(probability))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end
	
	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	elseif !accept && updateHist #this is the only change compared to HistoMetC.
		MetC.hist[binOld] += 1
	end

	return accept

end


#=============================================================================#
#================================HISTOAbsMetC=================================#
#=============================================================================#

mutable struct HISTOAbsMetC <: MetC
	kT::Float64
	w::Float64
	delta::Float64
	resetPeriod::Float64
	refCNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOAbsMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOAbsMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOAbsMetC, oldCluster::Cluster, newCluster::Cluster)

	if newCluster.energy < oldCluster.energy
		return true
	end

	updateHist = false
	binNew = nothing
	put!(MetC.io[2], "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)")
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			put!(MetC.io[2], "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)")
		end
		MetC.timeElapsed += 1

		hScore = 0
	else
		histSum = sum(MetC.hist) # used to normalise bar heights
		if histSum == 0 # avoid division by zero
			histSum = 1
		end

		simOld = getCNASimilarity(getCNA(oldCluster), MetC.refCNA)
		binOld = simOld == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simOld/MetC.delta) + 1 # get bin of old Cluster
		hOld   = MetC.hist[binOld]/histSum # get height (normalised) of bars
		
		put!(MetC.io[2], "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)")

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		put!(MetC.io[2], "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)")

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = (1 - MetC.w) * exp((oldCluster.energy - newCluster.energy) / MetC.kT) + hScore

	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(MetC.io[2], "\nhScore = $(hScore)\nChance to accept = " * string(probability))
	while isready(MetC.io[2])
		print(MetC.io[1], take!(MetC.io[2]))
	end
	
	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	end

	return accept

end
