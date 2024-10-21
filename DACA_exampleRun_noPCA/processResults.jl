include("/Users/smini250/Documents/DACGround/v1.2.4/DAC/src/DAC.jl")

import .DAC
const DC = DAC

using JLD2
using BenchmarkTools
using PyPlot

function getAtom80Class(clusterVector::DC.ClusterVector)
	# list of all atomic classes.
    allClasses = Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, 0)
    # list of class frequencies
    freq = Vector{Int64}(undef, 0)

    # for each structure
    for i in 1:clusterVector.N[]
        # calculate CNA
        nCNA::DC.normalCNAProfile = DC.getNormalCNAProfile(clusterVector.vec[i].positions, rcut)

        #for each atom in the cna profile
        for j in 1:length(nCNA)
            # if atomic class is new, add an entry for it
            if !(nCNA[j] in allClasses)
                push!(allClasses, nCNA[j])
                push!(freq, 0)
            # if atomic class already known, increment frequency
            else
                index = findfirst(x->x==nCNA[j], allClasses)
                freq[index] += 1
            end
        end
    end

    # get sorted order of frequencies from largest to smallest
    p = sortperm(freq, rev=true)
    # get total count of frequencies
    total = sum(freq)
    # get proportions of each class
    prop = [i/total for i in freq]

    # add most frequent atomic classes to `classes' if not already present (i.e. not already in Atom-64-Class).
    # do this until the desired number of classes is reached.
    i = 1
    while length(classes) < 80
        if !(allClasses[p[i]] in classes)
            push!(classes, allClasses[p[i]])
        end

        i += 1
    end

    # make the final class the `undefined' class (i.e. none of the above.)
    temp = classes[80]
    classes[80] = classes[64]
    classes[64] = temp
	
	return classes
end

function main()
	# read the clusterVector. filename is argument of jldopen function. name of datastructure is the 2nd string
	# which is always stored as "clusterVector"
	clusterVector = jldopen("clusterVector.jld2")["clusterVector"]
	rcut = 3.4765 # rcut for calculating CNA signatures / defining minimum bond length

	#get the number of structures stored in clusterVector
	N = clusterVector.N[]

	#the positions, energies, and CNA profiles of the structures can be accessed the following way.
	#here this is done for strucure 1, the lowest energy structure present in the clusterVector.
	#(structures are stored in ascending order of energies in this data-structure).
	clusterVector.vec[1].energy # energy of structure 1
	clusterVector.vec[1].positions # atomic coordinates of structure 1
	clusterVector.vec[1].CNA # CNA profile of structure 1

	# to retrieve the energies and similarity compared to a reference structure (in this
	# case the lowest energy structure), the following function can be used.
	# this data can then either be plotted with Julia (e.g. using the PyPlot package)
	# or alternatively printed for use in other software.
	# plotting sims on the x-axis and energies on the y-axis yields the 
	# energy vs similarity plots present in the Paper.
	sims, energies = DC.getSimsAndEnergies(clusterVector, clusterVector.vec[1].CNA)


	# the getClasses function can be used to retrieve the atom classes used in Atom-64-Class
	classes = DC.getClasses()
	
	# the below is an alternate is Atom-80-Class was used instead.
	#classes = getAtom80Class()

	# the getClassMatrix function can be used to calculate the frequency of a given set of atom
	# classes within a clusterVector.
	# the function returns a Matrix of size (<number of classes>, <number of structures>).
	# classMatrix[1, 1] then would give the count of class 1 in structure 1 of the clusterVector.
	classMatrix = DC.getClassMatrix(clusterVector, rcut, classes)


	# the following can be used to retrieve the gaussian mixture model and, if used,
	# principal component analysis model, that were trained when the algorithm was run.
	gmmpca = jldopen("GMM_PCA.jld2")
	gmm = gmmpca["gmm"]
	#pca = gmmpca["pca"]


	# the following method can be used to determine which cluster (classes) a set of structures belong to
	#structureClasses = DC.getStructureClasses(clusterVector, gmm, pca, rcut, classes)

	# if PCA was not used, the following method can be used instead to determine cluster assignments.
	classMatrixFloat::Matrix{Float64} = classMatrix'[:, :]
	structureClasses = DC.getStructureClasses(clusterVector, gmm, classMatrixFloat)
	
	# the below is used to classify NPs as ICO, PICO, AICO, TWI, FCC, DEC, or AMB.
	# note an additional class is used internally by our group (EXICO). In the paper,
	# EXICO is presented at AICO (ANTIICO).
	# HCP is also a possibility but no such structures were encountered in the paper.
	# TRICO and TWICO were classified as PICO (POLYICO) in the paper.
	nanoparticleClassifications = Vector{String}(undef, N)
	for i in 1:N
		nanoparticleClassifications = DC.classifyCluster(clusterVector.vec[i].positions, rcut)
	end

end

main()
