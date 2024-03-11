
"""
	getSimsAndEnergies(clusterVector::ClusterVector, refCNA::CNAProfile)

returns the energies and sims from a ClusterVector relative to a given 
reference CNA profile.
"""
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

"""
	getSimsAndEnergiesAndClassMatrix(clusterVector::ClusterVector, rcut::Float64)

returns the energies and sims from a ClusterVector relative to a given 
reference CNA profile, as well as the atom classes (frequencies) from trained GMM and PCA models. 
rcut is required to recalculate normal CNA profiles.
"""
function getClassMatrix(clusterVector::ClusterVector, rcut::Float64)
	nSamples::Int64 = clusterVector.N[]

	# class definitions
	classes = getClasses()
	nClasses::Int64 = length(classes)

	# class of each structure in clusterVector
	structureClasses = Vector{Int64}(undef, nSamples)
	# atomic classes
	atomClassMatrix = Matrix{UInt8}(undef, nClasses, nSamples)

	for i in 1:nSamples
		# get normal CNA, atomic classes then frequency of each atomic class
		nCNA::normalCNAProfile = getNormalCNAProfile(clusterVector.vec[i].positions, rcut)
		atomClasses = getAtomClasses(nCNA, classes)
		atomClassMatrix[:, i] = getFrequencyClassVector(atomClasses, nClasses, UInt8)
	end

	return atomClassMatrix
end


"""
	getSimsAndEnergiesAndPCAs(clusterVector::ClusterVector, pca::PCA, rcut::Float64)

returns the energies and sims from a ClusterVector relative to a given 
reference CNA profile, as well as the PCA data for a given PCA model, with a given rcut
to recalculate normal CNA profiles.
"""
function getPCAxes(clusterVector::ClusterVector, pca::PCA, atomClassMatrix::Matrix{UInt8})

	nSamples::Int64 = clusterVector.N[]

	# convert atom class frequency to PC values
	X = Matrix{Float64}(undef, nSamples, size(pca)[2])
	X[:, :] = predict(pca, atomClassMatrix)'[:, :]

	return X
end


getPCAxes(clusterVector::ClusterVector, pca::PCA, rcut::Float64) = getPCAxes(clusterVector, pca, getClassMatrix(clusterVector, rcut))


"""
	getStructureClasses(clusterVector::ClusterVector, gmm::GMM, X::Matrix{Float64})

returns the classes from trained GMM model given the PCA space data 
"""
function getStructureClasses(clusterVector::ClusterVector, gmm::GMM, X::Matrix{Float64})

	nSamples::Int64 = clusterVector.N[]

	# get the division index of each structure
	posteriorprobs = gmmposterior(gmm, X)[1]
	divisions = findmax(posteriorprobs, dims=2)[2]
	structureClasses = Vector{Int64}(undef, nSamples)
	for i in 1:nSamples
		structureClasses[i] = divisions[i][2]
	end

	return structureClasses

end

"""
	getStructureClasses(clusterVector::ClusterVector, gmm::GMM, pca::PCA, rcut::Float64)

returns the structure classes from a trained GMM model. Wrapper function for an 
input of the trained PCA model and rCut
"""
getStructureClasses(clusterVector::ClusterVector, gmm::GMM, pca::PCA, rcut::Float64) = getStructureClasses(clusterVector, gmm, getPCAxes(clusterVector, pca, rcut))


function getAxesLims(system::String)
	if system == "Au55"
		return [-195.4, -192.0]
	elseif system == "LJ75"
		return [-397.5, -360.0]
	elseif system == "LJ98"
		return [-543.7, -500.0]
	end
end


####################################################
####################################################
####################plotBirdpoo#####################
####################################################
####################################################

"""
	plotBirdpoo(sims::Vector{Float64}, energies::Vector{Float64}, system::String, c::Union{Vector{Int64}, Nothing}=nothing, filename::String="")

Takes a Vector of sims and energies and plots a birdpoo plot. system required to set ylim.
colours can be fed, e.g. from structure classes from GMM. 
if filename is provided, saves the figure, otherwise figure is displated with `show()`.
"""
function plotBirdpoo(sims::Vector{Float64}, energies::Vector{Float64}, system::String, c::Union{Vector{Int64}, Nothing}=nothing, filename::String="")
	
	display = filename == ""

	scatter(sims, energies, c=c, s=1)

	xlim([0, 1])

	ylim(getAxesLims(system))

	if display
		show()
	else
		savefig(filename, dpi=250)
	end

	return nothing
end

"""
	plotBirdpoo(clusterVector::ClusterVector, refCNA::CNAProfile, gmm::GMM, pca::PCA, rcut::Float64, system::String, filename::String="")  

Wrapper function for plotting a birdpoo plot. Takes a clusterVector, and refCNA to get sims and energies from. Also takes
gmm and pca models, alongside rcut, to recalculate normal cna profiles to determine structure classes from gmm model.
passes these values along side system and filename to the actual plotting function.
"""
function plotBirdpoo(clusterVector::ClusterVector, refCNA::CNAProfile, gmm::GMM, pca::PCA, rcut::Float64, system::String, filename::String="")  

	sims, energies = getSimsAndEnergies(clusterVector, refCNA)
	structureClasses = getStructureClasses(clusterVector, gmm, pca, rcut)

	plotBirdpoo(sims, energies, system, structureClasses, filename)

	return nothing
end

"""
	plotBirdpoo(clusterVector::String, refCNA::String, gmm::String, pca::String,
											rcut::Float64, system::String, filename::String="")

Wrapper function for plotting birdpoo plot. Takes filenames for clusterVector, structurename of refCNA, and 
filenames of gmm and pca models (all files names point to .jld2 file), opens the files and passes to 
middleman processing function to get plotting data.
"""
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
																													nothing,
																													filename)

plotBirdpoo(clusterVector::String, refCNA::CNAProfile, system::String, filename::String="") = plotBirdpoo(jldopen(clusterVector)["clusterVector"], 
																											refCNA, 
																											system, 
																											filename)

plotBirdpoo(clusterVector::String, refCNA::String, system::String, filename::String="") = plotBirdpoo(jldopen(clusterVector)["clusterVector"], 
																										stringToCNA(getCNA(refCNA)), 
																										system, 
																										filename)










####################################################
####################################################
#############plotBirdpooAndILSDistances#############
####################################################
####################################################


function plotBirdpooAndILSDistances(sims::Vector{Float64}, energies::Vector{Float64}, Ri::Vector{Float64}, iterationLabelledAt::Vector{Int64}, system::String; filename::String="", cmap::String="jet")

	display = filename == ""

	fig, axs = subplots(1, 2)

	x = [i for i in 1:length(Ri)]
	N = length(Ri)
	axs[1].scatter(x, Ri, c=x, cmap=cmap)
	axs[1].plot(x, Ri, c="k")

	insert!(iterationLabelledAt, 1, 0)
	axs[2].scatter(sims, energies, c=iterationLabelledAt, s=1, cmap=cmap)

	axs[2].set_xlim([0, 1])

	axs[2].set_ylim(getAxesLims(system))

	if display
		show()
	else
		savefig(filename, dpi=250)
	end

	return nothing
end


function plotBirdpooAndILSDistances(clusterVector::ClusterVector, refCNA::CNAProfile, rcut::Float64, 
									system::String; filename::String="", cmap::String="jet", cutOff::Int64=-1)
	if cutOff != -1
		clusterVector.vec = clusterVector.vec[1:cutOff]
		clusterVector.N = Threads.Atomic{Int64}(cutOff)
	end

	
	sims, energies = getSimsAndEnergies(clusterVector, refCNA)
	classMatrix = getClassMatrix(clusterVector, rcut)
	N = length(sims)
	labels = zeros(Int64, N)
	labels[1] = 1
	Ri, iterationLabelledAt = ILS(classMatrix, labels, true)

	plotBirdpooAndILSDistances(sims, energies, Ri, iterationLabelledAt, system, filename=filename)

	return nothing
end

plotBirdpooAndILSDistances(clusterVector::String, refCNA::String, rcut::Float64, system::String; 
							filename::String="", cmap::String="jet", cutOff::Int64=-1) = plotBirdpooAndILSDistances(jldopen(clusterVector)["clusterVector"],
																													stringToCNA(getCNA(refCNA)),
																													rcut,
																													system,
																													filename=filename,
																													cmap=cmap,
																													cutOff=cutOff)
	









####################################################
####################################################
##################plotBirdpooAndPCA#################
####################################################
####################################################

function isNum(a::String)
    return tryparse(Int64, a) !== nothing
end

function plotBirdpooAndPCA(sims::Vector{Float64}, energies::Vector{Float64}, pcAxes::Any,
							plotGridSpecs::Tuple{Int64, Int64}, plotAxes::Vector{Tuple{Any, Any}}, system::String; 
							c::Union{Vector{Int64}, Nothing}=nothing, cmap::String="tab20", filename::String="")

	display = filename == ""

	fig, axs = subplots(plotGridSpecs...)
	nPlots::Int64 = plotGridSpecs[1] * plotGridSpecs[2]
	if length(plotAxes) != nPlots
		error("The number of specified plots in `plotAxes` must match the product of plotGridSpecs.")
	end

	for i in 1:nPlots
		x = nothing
		y = nothing

		if typeof(plotAxes[i][1]) == Int64
			x = pcAxes[:, plotAxes[i][1]]
			axs[i].set_xlabel("PC$(plotAxes[i][1])")
		elseif plotAxes[i][1] == "e"
			x = energies
			axs[i].set_xlim(getAxesLims(system))
			axs[i].set_xlabel("Energy")
		elseif plotAxes[i][1] == "s"
			x = sims
			axs[i].set_xlim([0, 1])
			axs[i].set_xlabel("Similarity")
		end

		if typeof(plotAxes[i][2]) == Int64
			y = pcAxes[:, plotAxes[i][2]]
			axs[i].set_ylabel("PC$(plotAxes[i][2])")
		elseif plotAxes[i][2] == "e"
			y = energies
			axs[i].set_ylim(getAxesLims(system))
			axs[i].set_ylabel("Energy")
		elseif plotAxes[i][2] == "s"
			y = sims
			axs[i].set_ylim([0, 1])
			axs[i].set_ylabel("Similarity")
		end

		axs[i].scatter(x, y, s=1, c=c)

	end

	if display
		show()
	else
		savefig(filename, dpi=250)
	end

	return nothing

end

function plotBirdpooAndPCA(clusterVector::ClusterVector, refCNA::CNAProfile, pca::PCA, rcut::Float64,
							plotGridSpecs::Tuple{Int64, Int64}, plotAxes::Vector{Tuple{Any, Any}}, system::String; 
							gmm::Union{GMM, Nothing}=nothing, cmap::String="tab20", cutOff::Int64=-1, filename::String="")
	
	if cutOff != -1
		clusterVector.vec = clusterVector.vec[1:cutOff]
		clusterVector.N = Threads.Atomic{Int64}(cutOff)
	end

	sims, energies = getSimsAndEnergies(clusterVector, refCNA)
	pcAxes = getPCAxes(clusterVector, pca, rcut)
	if gmm != nothing
		structureClasses = getStructureClasses(clusterVector, gmm, pcAxes)
	end

	plotBirdpooAndPCA(sims, energies, pcAxes, plotGridSpecs, plotAxes, system, c=structureClasses, cmap=cmap, filename=filename)

	return nothing
end

plotBirdpooAndPCA(clusterVector::String, refCNA::String, pca::String, rcut::Float64, 
					plotGridSpecs::Tuple{Int64, Int64}, plotAxes::Vector{Tuple{Any, Any}}, system::String; 
					cmap::String="tab20", gmm::Union{String, Nothing}=nothing, cutOff::Int64=-1, filename::String="") = plotBirdpooAndPCA(jldopen(clusterVector)["clusterVector"],
																														stringToCNA(getCNA(refCNA)),
																														jldopen(pca)["pca"],
																														rcut,
																														plotGridSpecs,
																														plotAxes,
																														system,
																														cmap=cmap,
																														cutOff=cutOff,
																														gmm=gmm==nothing ? nothing : jldopen(gmm)["gmm"],
																														filename=filename)

