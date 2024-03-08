

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

function getSimsAndEnergiesAndClassMatrix(clusterVector::ClusterVector, refCNA::CNAProfile, rcut::Float64)
	nSamples::Int64 = clusterVector.N[]
	energies = Vector{Float64}(undef, nSamples)
	sims = Vector{Float64}(undef, nSamples)

	# class definitions
	classes = getClasses()
	nClasses::Int64 = length(classes)

	# class of each structure in clusterVector
	structureClasses = Vector{Int64}(undef, nSamples)
	# atomic classes
	atomClassMatrix = Matrix{Float64}(undef, nClasses, nSamples)

	for i in 1:nSamples
		#get energy and sim
		energies[i] = clusterVector.vec[i].energy
		sims[i] = getCNASimilarity(clusterVector.vec[i].CNA, refCNA)

		# get normal CNA, atomic classes then frequency of each atomic class
		nCNA::normalCNAProfile = getNormalCNAProfile(clusterVector.vec[i].positions, rcut)
		atomClasses = getAtomClasses(nCNA, classes)
		atomClassMatrix[:, i] = getFrequencyClassVector(atomClasses, nClasses, UInt8)
	end

	return sims, energies, atomClassMatrix
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


function plotBirdpoo(clusterVector::ClusterVector, refCNA::CNAProfile, gmm::GMM, pca::PCA, rcut::Float64, system::String, filename::String="")  
	sims, energies, structureClasses = getSimsAndEnergiesAndClasses(clusterVector, refCNA, gmm, pca, rcut)
	return plotBirdpoo(sims, energies, system, structureClasses, filename)
end


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


function plotBirdpooAndILSDistances(sims::Vector{Float64}, energies::Vector{Float64}, Ri::Vector{Float64}, iterationLabelledAt::Vector{Int64}, system::String, filename::String="")

	display = filename == ""

	fig, axs = subplots(1, 2)

	axs[1].scatter(sims, energies, c=iterationLabelledAt, s=1)

	axs[1].set_xlim([0, 1])

	if system == "Au55"
		axs[1].set_ylim([-195.4, -192.0])
	elseif system == "LJ75"
		axs[1].set_ylim([-397.5, -360.0])
	end

	x = [i for i in 1:length(Ri)]
	N = length(Ri)
	axs[2].scatter(x, Ri, c=iterationLabelledAt)
	axs[2].plot(x, Ri, c="k")

	if display
		show()
	else
		savefig(filename, dpi=250)
	end
end


function plotBirdpooAndILSDistances(clusterVector::String, refCNA::String, rcut::Float64, system::String, filename::String="", cutOff::Int64=-1)
	clusterVector = jldopen(clusterVector)["clusterVector"]
	refCNA = stringToCNA(getCNA(refCNA))
	if cutOff != -1
		clusterVector.vec = clusterVector.vec[1:cutOff]
	end

	
	sims, energies, classMatrix = getSimsAndEnergiesAndClassMatrix(clusterVector, refCNA, rcut)
	N = length(sims)
	labels = zeros(Int64, N)
	labels[1] = 1
	Ri, iterationLabelledAt = ILS(classMatrix, labels, true)

	return plotBirdpooAndILSDistances(sims, energies, Ri, iterationLabelledAt, system, filename)
end

