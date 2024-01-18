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
	useExplorationDataOnly::Bool
	io::Tuple{IO, Channel}
end

function getAcceptanceBoolean(MetC::ITCMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		accept = true
	else
		probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
		
		metcLog *= "\nChance to accept = $(string(probability))"
		
		accept = probability > rand()
	end

	# if the hop is rejected before any GMM checks are made, stop here
	if !accept
		return accept, metcLog
	end


	sim::Float64 = getCNASimilarity(MetC.refCNA, newCluster.CNA)
	metcLog *= "\nnewCluster is $sim similar to reference."
	if sim < MetC.threshold
		return false, metcLog
	end

	return accept, metcLog
end

#=============================================================================#
#==================================GMMMetC====================================#
#=============================================================================#

mutable struct GMMMetC <: MetC
	gaussian::GMM
	gaussianCluster::Int64
	pca::PCA
	mode::Symbol
	useExplorationDataOnly::Bool
	kT::Float64
	classes::normalCNAProfile
	nClasses::Int64
	workspace::Matrix{Float64}
	io::Tuple{IO, Channel}
end

function GMMMetC(gaussian::GMM, gaussianCluster::Int64, pca::PCA, mode::Symbol, useExplorationDataOnly::Bool, kT::Float64, io::Tuple{IO, Channel})
	# sets workspace as a 1x{PCA_out_dims} Matrix.
	classes = getClasses()
	GMMMetC(gaussian, gaussianCluster, pca, mode, useExplorationDataOnly, kT, classes, length(classes), Matrix{Float64}(undef, 1, size(pca)[2]), io)
end

function setMLClusterIndex!(MetC::GMMMetC, cluster::Cluster)
	fractionalClassVector = getFrequencyClassVector(getAtomClasses(cluster.nCNA, MetC.classes), MetC.nClasses)
	bh.metC.workspace[1, :] = predict(MetC.pca, fractionalClassVector)'[:, :]
	mlClusterIndex = findmax(gmmposterior(MetC.gaussian, MetC.workspace)[1])[2][2]
end

"""
	getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on the EnergyMetC.
"""
function getAcceptanceBoolean(MetC::GMMMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	if newCluster.energy < oldCluster.energy
		accept = true
	else

		probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
		
		metcLog *= "\nChance to accept = $(string(probability))"
		
		accept = probability > rand()
	end

	# if the hop is rejected before any GMM checks are made, stop here
	if !accept
		return accept, metcLog
	end

	# get the class vector for atom classes (Roncaglia scheme)
	fractionalClassVector = newCluster.atomClassCount

	# transform the class vector into PCA space
	MetC.workspace[1, :] = predict(MetC.pca, fractionalClassVector)'[:, :]

	# get the probabilities that this datapoint belongs to each of the n Gaussian clusters.
	posteriorProbs = gmmposterior(MetC.gaussian, MetC.workspace)[1]

	# this mode only accepts a hop if the target Gaussian cluster is the most likely Gaussian for this datapoint
	if MetC.mode == :maxProbOnly
		# `findmax` returns (maxvalue, indexOf), where indexOf is of type CartesianIndex{2} (as the arg is a 1xn Matrix).
		mlClusterIndex = findmax(posteriorProbs)[2][2]
		accept = mlClusterIndex == MetC.gaussianCluster
		metcLog *= "\nnewCluster belongs to cluster $(mlClusterIndex)."
		return accept, metcLog
	# this accepts a hop depending on how likely it is this datapoint belongs to the target Gaussian.
	elseif MetC.mode == :clusterProb
		accept = posteriorProbs[1, MetC.gaussianCluster] > rand()
		metcLog *= "\nnewCluster belongs to cluster $(MetC.gaussianCluster)\n\twith probability $(posteriorProbs[1, MetC.gaussianCluster])"
		return accept, metcLog
	end

	return accept, metcLog * "__badReturn"

end

#=============================================================================#
#=============================GMMwithInfTempMetC==============================#
#=============================================================================#

mutable struct GMMwithInfTempMetC <: MetC
	gaussian::GMM
	gaussianCluster::Int64
	pca::PCA
	mode::Symbol
	useExplorationDataOnly::Bool
	kT::Float64
	classes::normalCNAProfile
	nClasses::Int64
	workspace::Matrix{Float64}
	infTempPeriod::Int64
	infTempDuration::Int64
	hopsToInfTemp::Int64
	hopsRemainingOfInfTemp::Int64
	infTempEnergyToBeat::Float64
	io::Tuple{IO, Channel}
end

function GMMwithInfTempMetC(gaussian::GMM, gaussianCluster::Int64, pca::PCA, mode::Symbol, useExplorationDataOnly::Bool, kT::Float64, infTempPeriod::Int64, 
							infTempDuration::Int64, hopsToInfTemp::Int64, hopsRemainingOfInfTemp::Int64, infTempEnergyToBeat::Float64, io::Tuple{IO, Channel})
	# sets workspace as a 1x{PCA_out_dims} Matrix.
	classes = getClasses()
	GMMwithInfTempMetC(gaussian, gaussianCluster, pca, mode, useExplorationDataOnly, kT, classes, length(classes), Matrix{Float64}(undef, 1, size(pca)[2]), 
					   infTempPeriod, infTempDuration, hopsToInfTemp, hopsRemainingOfInfTemp, infTempEnergyToBeat, io)
end

function setMLClusterIndex!(MetC::GMMwithInfTempMetC, cluster::Cluster)
	fractionalClassVector = getFrequencyClassVector(getAtomClasses(cluster.nCNA, MetC.classes), MetC.nClasses)
	bh.metC.workspace[1, :] = predict(MetC.pca, fractionalClassVector)'[:, :]
	mlClusterIndex = findmax(gmmposterior(MetC.gaussian, MetC.workspace)[1])[2][2]
end

"""
	getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on the EnergyMetC.
"""
function getAcceptanceBoolean(MetC::GMMwithInfTempMetC, oldCluster::Cluster, newCluster::Cluster)
	metcLog = ""
	
	tempIsInf = MetC.hopsRemainingOfInfTemp > 0 # is the temperature infinite?
	if tempIsInf
		metcLog *= "\nTemperature is infinite. Auto accepting uphill moves except for GMM conditions."

	# if kT is currently infinite OR if move is downhill, accept hop
	if tempIsInf || newCluster.energy < oldCluster.energy
		accept = true
		
	else # apply boltzmann to get chance to accept uphill move 
		probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
		
		metcLog *= "\nChance to accept = $(string(probability))"
		
		accept = probability > rand()
	end



	# if the hop is rejected before any GMM checks are made, stop here
	if !accept
		return accept, metcLog
	end

	# get the class vector for atom classes (Roncaglia scheme)
	fractionalClassVector = newCluster.atomClassCount

	# transform the class vector into PCA space
	MetC.workspace[1, :] = predict(MetC.pca, fractionalClassVector)'[:, :]

	# get the probabilities that this datapoint belongs to each of the n Gaussian clusters.
	posteriorProbs = gmmposterior(MetC.gaussian, MetC.workspace)[1]

	# this mode only accepts a hop if the target Gaussian cluster is the most likely Gaussian for this datapoint
	if MetC.mode == :maxProbOnly
		# `findmax` returns (maxvalue, indexOf), where indexOf is of type CartesianIndex{2} (as the arg is a 1xn Matrix).
		mlClusterIndex = findmax(posteriorProbs)[2][2]
		accept = mlClusterIndex == MetC.gaussianCluster
		metcLog *= "\nnewCluster belongs to cluster $(mlClusterIndex)."
	# this accepts a hop depending on how likely it is this datapoint belongs to the target Gaussian.
	elseif MetC.mode == :clusterProb
		accept = posteriorProbs[1, MetC.gaussianCluster] > rand()
		metcLog *= "\nnewCluster belongs to cluster $(MetC.gaussianCluster)\n\twith probability $(posteriorProbs[1, MetC.gaussianCluster])"
	end

	# if temperature is NOT currently infinite
	if !tempIsInf 
		if accept && newCluster.energy < MetC.infTempEnergyToBeat # if a new energy minimum has been found and hopped to since the last inf temp run
			MetC.hopsToInfTemp = MetC.infTempPeriod # reset the hops until infinite temperature
		else # decrement remaining hops until infinite temperature
			MetC.hopsToInfTemp -= 1
			if MetC.hopsToInfTemp == 0 # if NEXT hop IS at infinite temperature, set the number of hops of infinite temperature.
				MetC.hopsRemainingOfInfTemp = MetC.infTempDuration
			end
		end

	# otherwise if temperature IS currently infinite, decrement remaining hops of infinite temperature
	else 
		MetC.hopsRemainingOfInfTemp -= 1
		if MetC.hopsRemainingOfInfTemp == 0 # if this was the last hop of infinite temperature...
			MetC.hopsToInfTemp = MetC.infTempPeriod # reset number of hops until next reseed period
			MetC.infTempEnergyToBeat = accept ? newCluster.energy : oldCluster.energy # set the energy to beat depending on if the hop was accepted.
		end

	end

	return accept, metcLog

end
