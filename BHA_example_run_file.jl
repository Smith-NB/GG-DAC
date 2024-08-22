include("<path_to_dir_of_DAC.jl>/DAC.jl")
import .DAC
const DC = DAC

using JLD2
using BenchmarkTools

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
    fmaxTight = 0.001 # tighter convergence criterion
    tightEnergyThreshold = -195.18 # threshold below which tighter convergence is used
    rcut = 3.4765 # rcut for calculating CNA signatures / defining minimum bond length
    energyRounding = 3 # dp that energies are stored in clusterVector with. 3 used for Au55, 2 for LJs.

    kT = 0.1            # temperature
    walltime = 3.5 # in hours
    recordingMode = "all"

    #output log files.
    io = (open("DAC.out", "w"), Channel{String}(1))
    logio = (open("log.txt", "w"), Channel{String}(1))
    cnaio = (open("CNAlog.txt", "w"), Channel{String}(1))


    # for storing structural data
    clusterVector = DC.ClusterVector(Vector{DC.ClusterCompressed}(), Threads.Atomic{Int64}(0), ReentrantLock())	

    # n_reseed (number of hops that can be performed since the last improvement to before a reseed is triggered)
    # 70, 100, 100, used for reseedPeriod for LJ75, 98 and 104.
    reseedPeriod = 50

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

    logResumeFile = true
    # do not start a new walk on a reseed for BHA.
    exitOnReseed = false
    # number of minimisations
    steps = 310000
    stepsAtomic = Threads.Atomic{Int64}(0)

    # seed generated randomly.
    seed = DC.generateRandomSeed(formula, boxLength, vacuumAdd, false, coherencyDistance)

    # what set of tasks to perform after each optimisation. A set of options is provided in BasinHopping.jl
    postOptimisationTasks = DC.standardPostOptimisationTasks!
    
    # Metropolis Criterion (from MetC.jl)
    MetC = DC.EnergyMetC(kT,io) # Used for BHA
    
    # HISTO MetC used for DACA exploration phase.
    #w = 100.0 # weight given to HISTO memory contribution
    #delta = 0.10 # histogram bar widths
    #waitTime = 100 # number of minimisations performed before reference set for SCM for HISTO feature space.
    #resetPeriod = Inf # redundant, leave as Inf.
    #MetC = DC.HISTOMetC(kT, w, delta, resetPeriod, clusterVector, waitTime, io)

    #reseeder
    Reseeder = DC.NewLESReseeder(reseedPeriod, reseedPeriod, Inf, DC.generateRandomSeed, [formula, boxLength, vacuumAdd, true, coherencyDistance])

    #energy potential
    rgl = DC.RGL(0.2197, 10.53, 4.30, 2.88, 1.855, N)

    # optimiser
    opt = DC.FIRE()
    wID = 1
    
    # clean up NP
    DC.setCalculator!(seed, rgl)
    DC.optimize!(opt, seed, fmax)

    # generate new structure if NP is not coherent (i.e. structure is 2+ seperate NPs.)
    while !DC.isClusterCoherent(seed.positions, coherencyDistance)
        DC.write_xyz("badCluster.xyz", seed)
        exit()
        seed = DC.generateRandomSeed(formula, boxLength, vacuumAdd, false, coherencyDistance)
        
        DC.centreCluster!(seed)
        DC.setCalculator!(seed, rgl)
        DC.optimize!(opt, seed, fmax)
    end

    ## create the walker (BasinHopping.jl).
    bh = DC.BasinHopper(opt, rgl, MetC, Reseeder, formula, boxLength, vacuumAdd, kT, perturber, postOptimisationTasks, fmax, fmaxTight,
                            tightEnergyThreshold, rcut, energyRounding, walltime, recordingMode, io, logio, cnaio, clusterVector, logResumeFile, 
                            exitOnReseed, "v1.2.4")
    # start the walk and store the number if minimisations performed.
    hopsPerformed = DC.hop(bh, steps, stepsAtomic, seed, wID, additionalInfo, start, "v1.2.4")

    # save structure data in JLD2 format.
    @save "clusterVector.jld2" clusterVector
    println("Size of clusterVector ", Base.summarysize(clusterVector), " $(clusterVector.N[])")
    close.(io)
    close.(logio)
    close.(cnaio)
end

main()
