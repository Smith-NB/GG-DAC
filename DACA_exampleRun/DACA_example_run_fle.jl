include("/Users/smini250/Documents/DACGround/v1.2.4/DAC/src/DAC.jl")
import .DAC
const DC = DAC

using JLD2
using BenchmarkTools

# checks if all clusters have been run for their allocated times. returns true if yes, else false.
function allClustersExploited(hopsPerCluster::Vector{Threads.Atomic{Int64}}, hopsToRun::Int64)
    for i in 1:length(hopsPerCluster)
        if hopsPerCluster[i][] < hopsToRun
            return false
        end
    end

    return true
end

function main()
    start = DC.now() # get start time
    N = 55 # number of atoms
    formula = Dict("Au" => N) # nanoparticle formula. "Ne" used for LJ systems.
    formulaString = "" 
    for k in keys(formula)
        formulaString *= "$k"
        n = formula[k]
        formulaString *= "$n"
    end

    boxLength = 5.0 # length of box to add atoms to for random NP structure generation
    vacuumAdd = 10.0 # vacuum to add to the box after structure generation
    fmax = 0.01     # optimisation convergence criterion

	# for Au55, the following setting is nessecary to ensure that the 16 lowest energy structures can be distinguished
	# based on their energies. 
	# To disable, set `fmaxTight` to `fmax`, and/or set `tightEnergyThreshold` to` -Inf`
    fmaxTight = 0.001 # tighter convergence criterion
    tightEnergyThreshold = -195.18 # threshold below which tighter convergence is used
    
	rcut = 3.4765 # rcut for calculating CNA signatures / defining minimum bond length
    energyRounding = 3 # dp that energies are stored in clusterVector with. 3 used for Au55, 2 for LJs.

    kT = 0.1            # temperature
    walltime = 3.5 # in hours
    recordingMode = "all"

    # output log files
    io = (open("DAC.out", "w"), Channel{String}(1))
    logio = (open("log.txt", "w"), Channel{String}(1))
    cnaio = (open("CNAlog.txt", "w"), Channel{String}(1))

    # for storing structural data
    clusterVector = DC.ClusterVector(Vector{DC.ClusterCompressed}(), Threads.Atomic{Int64}(0), ReentrantLock())	

    # n_reseed (number of hops that can be performed since the last improvement to before a reseed is triggered)
    # 500 used for Au55, 1000 for LJs
    reseedPeriod = 500 

    # atomic coordinates perturber
    function perturber(coords::Matrix{Float64})
        r = rand()
        if r < 0.75 # this is the likelihood of the CDO being used for perturbation
            return DC.perturbCluster(coords, 1.0), "CDO perturbation" # 1.0 is the stepwidth for the CDO (0.40 for LJs)
        end
        
        #print(io[1], "\nSAO perturbation")
        return DC.perturbClusterSurface(coords, 1, rcut), "SAO perturbation" # 1 is the number of atoms to move with SAO.
    end



    # for checking randomly generated structures are coherent (i.e. a single NP)
    coherencyDistance = 4.0 #4.0 for Aus, 2.0 for LJs
    additionalInfo::Dict{String, Any} = Dict()
    additionalInfo["coherencyDistance"] = coherencyDistance

    try 
        f = open("informationForResuming.txt", "r")
        lines = readlines(f)
        additionalInfo["stepsCompleted"] 		= parse(Int64,   split(lines[1], ":")[2])
        additionalInfo["Emin"] 					= parse(Float64, split(lines[2], ":")[2])
        additionalInfo["EminLocatedAt"] 		= parse(Int64,   split(lines[3], ":")[2])
        additionalInfo["reseedEnergyToBeat"] 	= parse(Float64, split(lines[4], ":")[2])
        additionalInfo["hopsToReseed"] 			= parse(Int64,   split(lines[5], ":")[2])
    catch SystemError 
        print("no resume file.\n") 	
    end

    logResumeFile = false # log resumption information.
    walkID = Threads.Atomic{Int64}(1) # tracks the number of walks performed
    nClusters = 4 # `k' for the GMM, number of clusters i.e. divisions.

    # the number of minimisations (aka hops) performed in each division/cluster
    clusterExploitationTimes = Vector{Threads.Atomic{Int64}}(undef, nClusters)
    for i in 1:nClusters 
        clusterExploitationTimes[i] = Threads.Atomic{Int64}(0) 
    end

    hopsToDoPerCluster = 100 #hops to perform per division. In this study, nClusters * hopsToDoPerCluster = 300000
    exitOnReseed = true # start a new walk after a reseed.

    # filename for exploration/training data.
    explorationDataFilename::String = ""
    for item in readdir()
        if startswith(item, "explorationData")
            explorationDataFilename = item
        end
    end
    explorationData = DC.jldopen(explorationDataFilename)["clusterVector"]

    # atomic classes defined by CNA signatures used as the feature space for training the gaussian mixture model
    classes = DC.getClasses()

    # the following commented section expands the above 64 classes to 80 classes by adding the 16 most prevalent classes from the training data
    # not already present in Atom-64-Class.
    # Uncomment this section to use Atom-80-Class instead of Atom-64-Class.
    #=

    # list of all atomic classes.
    allClasses = Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, 0)
    # list of class frequencies
    freq = Vector{Int64}(undef, 0)

    # for each structure
    for i in 1:explorationData.N[]
        # calculate CNA
        nCNA::DC.normalCNAProfile = DC.getNormalCNAProfile(explorationData.vec[i].positions, rcut)

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
    
    =#

    nClasses = length(classes) # number of classes.

    eLim_denominator = 4 # rho value for calculating E_thresh (Equation 3 in paper).

    # E_Thresh
    eLim = round(explorationData.vec[1].energy + (explorationData.vec[explorationData.N[]].energy - explorationData.vec[1].energy) / eLim_denominator, digits=2)

    nSamplesReduced::Int64 = 0 # number of training data structures with energy less than or equal to E_Thresh
    nSamplesTotal::Int64 = explorationData.N[] # total number of training structures
    for i in 1:explorationData.N[]
        if explorationData.vec[i].energy <= eLim
            nSamplesReduced += 1
        end
    end

    println("_________NSAMPLESREDUCED_____________ $nSamplesReduced")
    println("_________NSAMPLESTOTAL_____________ $nSamplesTotal")

    # matrix of percentage occurences of atomic classes for each structure in the training data
    classMatrixReduced = Matrix{Float64}(undef, nClasses, nSamplesReduced) # matrix for only structures with energy <= E_thresh. size = n_dims, n_samples
    classMatrixTotal = Matrix{Float64}(undef, nClasses, nSamplesTotal) # matrix for all training structures. size = n_dims, n_samples

    # for each training data structure...
    for i in 1:explorationData.N[]
        nCNA::DC.normalCNAProfile = DC.getNormalCNAProfile(explorationData.vec[i].positions, rcut) #calculate normal cna profile (cna signatures of each atom)
        atomClasses = DC.getAtomClasses(nCNA, classes) # get the classes of each atom
        classMatrixTotal[:, i] = DC.getFrequencyClassVector(atomClasses, nClasses, UInt8) # convert atomClasses to percentage occurences, and store in matrix
        # the line below is nessecary when PCA is disabled to prevent cholesky decomposition from failing in the Gaussian mixture model training.
        # It adds a small amount of random noise to each datapoint.
        #classMatrixTotal[:, i] += rand(nClasses)/100 
        if explorationData.vec[i].energy <= eLim # if energy below E_thresh, store in reduced matrix too.
            classMatrixReduced[:, i] = classMatrixTotal[:, i]
        end
    end
    
    # Use the following 4 lines if PCA is enabled
    pca = DC.fit(DC.PCA, classMatrixReduced; pratio=0.90) # train PCA model
    XReduced = DC.predict(pca, classMatrixReduced)'[:, :] # PCA transformed data (reduced)
    XTotal = DC.predict(pca, classMatrixTotal)'[:, :] # PCA transformed data (full)
    gmm = DC.GMM(nClusters, XReduced; kind=:full) # train GMM on reduced data.

    # Use the following line if PCA is disabled. The "'[:, :]" is to transpose the data.
    #gmm = DC.GMM(nClusters, classMatrixReduced'[:, :]; kind=:full)

    # probability of every reduced structure beloning to each cluster of GMM
    probsReduced = DC.gmmposterior(gmm, XReduced)[1] #PCA Enabled
    #probsReduced = DC.gmmposterior(gmm, classMatrixReduced'[:, :])[1] #PCA disabled

    # structures assigned to cluster they most likely belong to
    labelsReduced = [findmax(probsReduced[i, :])[2] for i in 1:nSamplesReduced]
    # calculate lowest energy structure of each cluster
    clusterLESs = [findfirst(i->i==j, labelsReduced) for j in 1:nClusters]
    # for each cluster, get the structures it contains
    clusterIndicesReduced = [findall(i->i==j, labelsReduced) for j in 1:nSamplesReduced]
    


    # probability of every structure beloning to each cluster of GMM
    probsTotal = DC.gmmposterior(gmm, XTotal)[1] #PCA Enabled
    #probsTotal = DC.gmmposterior(gmm, classMatrixTotal'[:, :])[1] #PCA Disabled
    # structures assigned to cluster they most likely belong to
    
    labelsTotal = [findmax(probsTotal[i, :])[2] for i in 1:nSamplesTotal]
    # for each cluster, get the structures it contains
    clusterIndicesTotal = [findall(i->i==j, labelsTotal) for j in 1:nSamplesTotal]

    nSeedSelections = 5 # upon a reseed, how many structures are drawn to choose the next walk's seed from (highest energy of these is selected).

    # stores the seed for the next walk of each cluster.
    seeds = Vector{DC.Cluster}(undef, nClusters)
    cell = zeros(Float64, 3, 3) # cell of the NP
    cell[1] = 20.0; cell[5] = 20.0; cell[9] = 20.0

    # set the seeds for the first walk in each cluster to the lowest energy structure of each cluster.
    for i in 1:nClusters
        seeds[i] = DC.Cluster(formula, explorationData.vec[clusterLESs[i]].positions, deepcopy(cell))
        DC.centreCluster!(seeds[i])
    end


    fillerTask = Task(begin end)
    schedule(fillerTask)
    tasks = fill(fillerTask, nClusters) # creates a Vector of Tasks that will return true for both istaskstarted and istaskdone

    tasklog = open("tasklog.txt", "w")

    nwhile = 0 # counter

    wIDs = Vector{Int64}(undef, nClusters) # walkIDs for tasklog

    # what set of tasks to perform after each optimisation. A set of options is provided in BasinHopping.jl
    postOptimisationTasks = DC.extendedPostOptimisationTasks! #PCA Enabled.
    #postOptimisationTasks = DC.extendedPostOptimisationTasksNoPCA! #PCA Disabled.

    # run until each cluster has had `hopsToDoPerCluster' minimisations performed.
    while !allClustersExploited(clusterExploitationTimes, hopsToDoPerCluster)
        # if walltime reached, break loop/
        if (DC.now() - start) / DC.Hour(1) > walltime
            break
        end
        
        # debugging aids
        println(tasklog, "$nwhile $((DC.now() - start) / DC.Hour(1) < walltime)\n"); 
        println(tasklog, tasks); 
        println(tasklog, wIDs); 
        for i in 1:nClusters print(tasklog, "$(clusterExploitationTimes[i][]) ") end
        print(tasklog, "\n")

        # for each cluster
        for k in 1:nClusters
            # if the previous walk is complete and there are still minimisations to perform in this cluster.
            if istaskdone(tasks[k]) && clusterExploitationTimes[k][] < hopsToDoPerCluster
                # spawn a new task (walk)
                tasks[k] = Threads.@spawn begin
                    
                    wID = Threads.atomic_add!(walkID, 1)
                    wIDs[k] = wID
                    println(tasklog, "starting task $wID on cluster $k")
                    
                    # Metropolis Criterion (from MetC.jl)
                    MetC = DC.GMMMetC(gmm, k, pca, :maxProbOnly, true, kT, classes, io) #PCA Enabled.
                    #MetC = DC.GMMnoPCAMetC(gmm, k, :maxProbOnly, true, kT, classes, ones(length(classes)), io) #PCA Disabled.
                    # Reseeder 
                    Reseeder = DC.NewLESReseeder(reseedPeriod, reseedPeriod, Inf, DC.generateRandomSeed, [formula, boxLength, vacuumAdd])

                    # Energy potential
                    rgl = DC.RGL(0.2197, 10.53, 4.30, 2.88, 1.855, N)

                    # optimiser
                    opt = DC.FIRE()
                    
                    # set the calculator of the seed.
                    DC.setCalculator!(seeds[k], rgl)
                    # optimize the seed.
                    DC.optimize!(opt, seeds[k], fmax)

                    # if the seed is incoherent, generate a new seed and make sure it is coherent.
                    while !DC.isClusterCoherent(seeds[k].positions, 4.0)
                        
                        # select `nSeedSelections' seeds and take the highest energy one.
                        seedIndex = rand(clusterIndicesTotal[k])
                        for j in 2:nSeedSelections
                            newSeedIndex = rand(clusterIndicesTotal[k])
                            if explorationData.vec[newSeedIndex].energy > explorationData.vec[seedIndex].energy
                                seedIndex = newSeedIndex
                            end
                        end
                        
                        # store the selected seed.
                        seeds[k] = DC.Cluster(formula, explorationData.vec[seedIndex].positions, cell)
                        
                        # clean up the seed.
                        DC.centreCluster!(seeds[k])
                        DC.setCalculator!(seeds[k], rgl)
                        DC.optimize!(opt, seeds[k], fmax)
                    end

                    # create the walker (BasinHopping.jl).
                    bh = DC.BasinHopper(opt, rgl, MetC, Reseeder, formula, boxLength, vacuumAdd, kT, perturber, postOptimisationTasks, fmax,
                                        fmaxTight, tightEnergyThreshold, rcut, energyRounding, walltime, recordingMode, io, logio, cnaio, 
                                        clusterVector, logResumeFile, exitOnReseed, "v1.2.4")

                    # start the walk and store the number if minimisations performed.
                    hopsPerformed = DC.hop(bh, hopsToDoPerCluster-clusterExploitationTimes[k][], clusterExploitationTimes[k], seeds[k], wID, additionalInfo, start, "v1.2.4")
                    
                    
                    println(tasklog, "completing task $wID on cluster $k")
                    
                     # select `nSeedSelections' seeds and take the highest energy one.
                    seedIndex = rand(clusterIndicesTotal[k])
                    for j in 2:nSeedSelections
                        newSeedIndex = rand(clusterIndicesTotal[k])
                        if explorationData.vec[newSeedIndex].energy > explorationData.vec[seedIndex].energy
                            seedIndex = newSeedIndex
                        end
                    end
                    # store the selected seed.
                    seeds[k] = DC.Cluster(formula, explorationData.vec[seedIndex].positions, cell)

                    # clean up the seed.
                    DC.centreCluster!(seeds[k])
                    println(tasklog, "ending task $wID on cluster $k")
                end
            end 
        end

        sleep(1)
    end

    while true
        tasksDone = true
        for k in 1:nClusters
            print("$(istaskdone(tasks[k])) ")
            if !istaskdone(tasks[k])
                tasksDone = false
                println("task $k incomplete: started $(istaskstarted(tasks[k])) completed $(istaskdone(tasks[k]))"); flush(stdout)
            end
        end
        println("\n$tasks"); flush(stdout)
        flush(io[1])
        flush(cnaio[1])
        flush(logio[1])
        if tasksDone
            break
        end
        sleep(1)
    end

    # save structure data in JLD2 format.
    @save "clusterVector.jld2" clusterVector
    # save the pca and gaussian mixture models in JLD2 format.
    jldsave("GMM_PCA.jld2"; gmm, pca)

    println("Size of clusterVector ", Base.summarysize(clusterVector), " $(clusterVector.N[])")
    close.(io)
    close.(logio)
    close.(cnaio)

end

#run
main()
