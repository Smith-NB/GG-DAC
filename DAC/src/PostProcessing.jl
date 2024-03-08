

function getSimsAndEnergies(clusterVector::ClusterVector, refCNA::CNAProfile)
	
	N = clusterVector.N[]
	energies = Vector{Float64}(undef, N)
	sims = Vector{Float64}(undef, N)
	for i in 1:N
		energies[i] = clusterVector.vec[i].energy
		sims[i] = getCNASimilarity(clusterVector.vec[i].CNA, refCNA)
	end

	return sims, energies

end

function getSimsAndEnergiesAndClasses(clusterVector::ClusterVector, refCNA::CNAProfile, gmm::GMM, pca::PCA, rcut::Float64)

	nSamples::Int64 = clusterVector.N[]
	energies = Vector{Float64}(undef, nSamples)
	sims = Vector{Float64}(undef, nSamples)

	# class definitions
	classes = getClasses()
	nClasses::Int64 = length(classes)

	# class of each structure in clusterVector
	structureClasses = Vector{Int64}(undef, nSamples)
	# atomic classes
	atomClassMatrix = Matrix{UInt8}(undef, nClasses, nSamples)

	for i in 1:nSamples
		#get energy and sim
		energies[i] = clusterVector.vec[i].energy
		sims[i] = getCNASimilarity(clusterVector.vec[i].CNA, refCNA)

		# get normal CNA, atomic classes then frequency of each atomic class
		nCNA::normalCNAProfile = getNormalCNAProfile(clusterVector.vec[i].positions, rcut)
		atomClasses = getAtomClasses(nCNA, classes)
		atomClassMatrix[:, i] = getFrequencyClassVector(atomClasses, nClasses, UInt8)
	end

	# convert atom class frequency to PC values
	X = Matrix{Float64}(undef, nSamples, size(pca)[2])
	X[:, :] = predict(pca, atomClassMatrix)'[:, :]

	# get the division index of each structure
	posteriorprobs = gmmposterior(gmm, X)[1]
	divisions = findmax(posteriorprobs, dims=2)[2]
	for i in 1:nSamples
		structureClasses[i] = divisions[i][2]
	end

	return sims, energies, structureClasses

end


function plotBirdpoo(sims::Vector{Float64}, energies::Vector{Float64}, system::String, c::Union{Vector{Int64}, Nothing}=nothing, filename::String="")
	
	display = filename == ""

	scatter(sims, energies, c=c, s=1)

	xlim([0, 1])

	if system == "Au55"
		ylim([-195.4, -192.0])
	elseif system == "LJ75"
		ylim([-397.5, -360.0])
	end

	if display
		show()
	else
		savefig(filename, dpi=250)
	end

end


plotBirdpoo(clusterVector::ClusterVector, refCNA::CNAProfile, gmm::GMM, pca::PCA, 
											rcut::Float64, system::String, filename::String="") = plotBirdpoo(getSimsAndEnergiesAndClasses(clusterVector, refCNA, gmm, pca, rcut)...,
																												system, filename)

plotBirdpoo(clusterVector::String, refCNA::String, gmm::String, pca::String,
											rcut::Float64, system::String, filename::String="") = plotBirdpoo(jldopen(clusterVector)["clusterVector"],
																												stringToCNA(getCNA(refCNA)),
																												jldopen(gmm)["gmm"],
																												jldopen(pca)["pca"],
																												rcut,
																												system,
																												filename)


plotBirdpoo(clusterVector::ClusterVector, refCNA::CNAProfile, system::String, filename::String="") = plotBirdpoo(getSimsAndEnergies(clusterVector, refCNA)...,  
																													system, 
																													filename)

plotBirdpoo(clusterVector::String, refCNA::CNAProfile, system::String, filename::String="") = plotBirdpoo(jldopen(clusterVector)["clusterVector"], 
																											refCNA, 
																											system, 
																											filename)

plotBirdpoo(clusterVector::String, refCNA::String, system::String, filename::String="") = plotBirdpoo(jldopen(clusterVector)["clusterVector"], 
																										stringToCNA(getCNA(refCNA)), 
																										system, 
																										filename)
