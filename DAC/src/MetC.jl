abstract type MetC end

struct EnergyMetC <: MetC
	kT::Float64
	io::Tuple{IO, Channel}
end

struct EnergyAndStructureMetC <: MetC
	cE::Float64
	cSCM::Float64
	eFunction::Function
	simFunction::Function
	autoAcceptDownhillMove::Bool
	io::Tuple{IO, Channel}
end

"""
	getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on  the EnergyMetC.
"""
function getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)
	
	if newCluster.energy < oldCluster.energy
		return true
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	
	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(io[2], "\nChance to accept = " * string(probability))
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
	
	
	accept = probability > rand()
	return accept

end


"""
	getAcceptanceBoolean(MetC::EnergyAndStructureMetC, oldCluster::Cluster, newCluster::Cluster)

Returns true or false for accepting the move from the oldCluster to the newCluster
	based on  the EnergyAndStructureMetC.
"""
function getAcceptanceBoolean(MetC::EnergyAndStructureMetC, oldCluster::Cluster, newCluster::Cluster)

	if MetC.autoAcceptDownhillMove && getEnergy(newCluster) < getEnergy(oldCluster)
		return true
	end

	eScore = MetC.eFunction(getEnergy(oldCluster), getEnergy(newCluster)) * MetC.cE
	simScore = MetC.simFunction(getCNASimilarity(getCNA(oldCluster), getCNA(newCluster))) * MetC.cSCM

	put!(io[2], "\nChance to accept = " * string(eScore + simScore))
	while isready(io[2])
		print(io[1], take!(io[2]))
	end

	accept = (eScore + simScore) > rand()
	return accept

end	