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
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		return true, metcLog
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	
	metcLog *= "\nChance to accept = $(string(probability))"
	
	
	accept = probability > rand()
	return accept, metcLog

end


#=============================================================================#
#======================EnergyWithReferenceRestrictionMetC=====================#
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
	metcLog = ""
	s = getCNASimilarity(MetC.refCNA, getCNAProfile(newCluster))

	if s < MetC.sigmaCut

		metcLog *=  "\nChance to accept = 0\nSimilarity to ref: $s"
		return false, metcLog
	end

	if newCluster.energy < oldCluster.energy
		return true, metcLog
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	
	# push a string to the channel then take the string in the channel and print it to the output file.
	metcLog *= "\nChance to accept = $(string(probability))"

	
	accept = probability > rand()
	return accept, metcLog

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
	metcLog = ""
	if MetC.autoAcceptDownhillMove && getEnergy(newCluster) < getEnergy(oldCluster)
		return true, metcLog
	end

	eScore = MetC.eFunction(getEnergy(oldCluster), getEnergy(newCluster)) * MetC.cE
	simScore = MetC.simFunction(getCNASimilarity(getCNA(oldCluster), getCNA(newCluster))) * MetC.cSCM

	metcLog *= "\nChance to accept = $(string(eScore + simScore))"


	accept = (eScore + simScore) > rand()
	return accept, metcLog

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
	refCNA::CNAProfile
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = CNAProfile()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		return true, metcLog
	end

	updateHist = false
	binNew = nothing
	metcLog *=  "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)"
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			metcLog *= "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)"
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
		
		metcLog *= "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)"

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		metcLog *= "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)"

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = exp((oldCluster.energy - newCluster.energy + hScore) / MetC.kT)	

	metcLog *=  "\nhScore = $(hScore)\nChance to accept = $(string(probability))"

	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	end

	return accept, metcLog

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
	refCNA::CNAProfile
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOAttMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = CNAProfile()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOAttMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOAttMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		return true, metcLog
	end

	updateHist = false
	binNew = nothing
	metcLog *= "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)"
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			metcLog *= "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)"
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
		
		metcLog *= "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)"

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		metcLog *= "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)"

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = exp((oldCluster.energy - newCluster.energy + hScore) / MetC.kT)	

	# push a string to the channel then take the string in the channel and print it to the output file.
	metcLog *= "\nhScore = $(hScore)\nChance to accept = $(string(probability))"

	
	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	elseif !accept && updateHist #this is the only change compared to HistoMetC.
		MetC.hist[binOld] += 1
	end

	return accept, metcLog

end


#=============================================================================#
#================================HISTOAbsMetC=================================#
#=============================================================================#

mutable struct HISTOAbsMetC <: MetC
	kT::Float64
	w::Float64
	delta::Float64
	resetPeriod::Float64
	refCNA::CNAProfile
	refID::Int64
	clusterVector::ClusterVector
	waitTime::Int64
	timeElapsed::Int64
	hist::Vector{Int64}
	io::Tuple{IO, Channel}
end

function HISTOAbsMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
	refCNA = CNAProfile()
	refID = -1
	timeElapsed = 0
	hist = zeros(Int64, trunc(Int64, 1/delta))
	return HISTOAbsMetC(kT, w, delta, resetPeriod, refCNA, refID, clusterVector, waitTime, timeElapsed, hist, io)
end

function getAcceptanceBoolean(MetC::HISTOAbsMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		return true, metcLog
	end

	updateHist = false
	binNew = nothing
	metcLog *= "\nlenRefCNA=$(length(MetC.refCNA)) timeElapsed=$(MetC.timeElapsed)"
	if length(MetC.refCNA) == 0
		if MetC.timeElapsed >= MetC.waitTime
			MetC.refCNA = MetC.clusterVector.vec[1].CNA
			MetC.refID = MetC.clusterVector.vec[1].ID
			metcLog *= "\nSETTING refCNA: $(MetC.refCNA)\nrefID: $(MetC.refID)"
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
		
		metcLog *= "\nsimOld=$(simOld)\nbinOld=$(binOld)\nhOld=$(hOld)"

		simNew = getCNASimilarity(getCNA(newCluster), MetC.refCNA)
		binNew = simNew == 1.0 ? trunc(Int64, 1/MetC.delta) : trunc(Int64, simNew/MetC.delta) + 1 # get bin of new Cluster
		hNew   = MetC.hist[binNew]/histSum # get height (normalised) of bars

		metcLog *= "\nsimNew=$(simNew)\nbinNew=$(binNew)\nhNew=$(hNew)"

		hScore = MetC.w * (hOld - hNew)
		updateHist = true
	end
	
	probability = (1 - MetC.w) * exp((oldCluster.energy - newCluster.energy) / MetC.kT) + hScore

	# push a string to the channel then take the string in the channel and print it to the output file.
	metcLog *= "\nhScore = $(hScore)\nChance to accept = $(string(probability))"
	
	accept = probability > rand()

	if accept && updateHist
		MetC.hist[binNew] += 1
	end

	return accept, metcLog

end

#=============================================================================#
#==================================ITCMetC====================================#
#=============================================================================#

mutable struct ITCMetC <: MetC
	kT::Float64
	threshold::Float64
	refCNA::CNAProfile
	explorationVector::ClusterVector
	useExplorationDataOnly::Bool
	nThreads::Int64
	io::Tuple{IO, Channel}
end

function getAcceptanceBoolean(MetC::ITCMetC, oldCluster::Cluster, newCluster::Cluster)

end

#=============================================================================#
#==================================GMMMetC====================================#
#=============================================================================#

mutable struct GMMMetC <: MetC
	kT::Float64
	threshold::Float64
	Gaussian::Gaussian
	GaussianCluster::Int64
	useExplorationDataOnly::Bool
	nThreads::Int64
	io::Tuple{IO, Channel}
end

function getAcceptanceBoolean(MetC::GMMMetC, oldCluster::Cluster, newCluster::Cluster)

end
