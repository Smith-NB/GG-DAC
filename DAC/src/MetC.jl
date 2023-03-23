abstract type MetC end

struct EnergyMetC <: MetC
	kT::Float64
	io::Tuple{IO, Channel}
end

function getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)
	
	if newCluster.energy < oldCluster.energy
		return true
	end

	probability = exp((oldCluster.energy - newCluster.energy) / MetC.kT)
	
	#=
	# push a string to the channel then take the string in the channel and print it to the output file.
	put!(io[2], "\nChance to accept = " * string(probability))
	while isready(io[2])
		print(io[1], take!(io[2]))
	end
	=#
	
	accept = probability > rand()
	return accept

end