#!/usr/bin/env julia

if "--help" in ARGS || "-h" in ARGS
    println()
    println("Usage: plot_birdpoo.jl [options]")
    println("Options:")
    println("  --fname <filename>   JLD2 file containing clusterVector (default: clusterVector.jld2)")
    println("  --ref <filename>     XYZ file to use as reference structure for CNA (default: first structure (LES) in clusterVector)")
    println("  --ref2 <filename>    XYZ file to use as second reference structure for CNA (default: off)")
    println("  --rcut <value>       Cutoff distance for CNA (default: 1.0). Can also be specified in BHA_run_script.jl or DACA_run_file.jl")
    println("                       only used if --ref is specified.")
    println("  -s, --schebarchov    Color points by Schebarchov class. Requires schebarchovClass.txt file in working directory")
    println("  -g, --gmm            Color points by GMM class. Requires gmmclasses.txt file in working directory")
    println("  -h, --help           Show this help message and exit")
    println()
    exit(0)
end

using DAC
import .DAC
const DC = DAC

using PyPlot

function get_rcut()

    rcut = 1.0
    runfilename = nothing
    if isfile("BHA_run_script.jl")
        runfilename = "BHA_run_script.jl"
    elseif isfile("DACA_run_file.jl")
        runfilename = "DACA_run_file.jl"
    end

    if "--rcut" in ARGS || "-r" in ARGS
        idx = findfirst(x -> x == "--rcut" || x == "-r", ARGS)
        if idx == length(ARGS)
            println("Error: --rcut requires an argument")
            exit(1)
        end
        rcut = parse(Float64, ARGS[idx+1])
        println("rcut = $rcut from command line arg")
    elseif runfilename != nothing
        lines = readlines(open(runfilename))
        global linei = 1
        while !occursin("rcut = ", lines[linei]) && linei < length(lines)
            global linei += 1
        end
        rcut = parse(Float64, split(lines[linei])[3])
        println("rcut = $rcut found in $runfilename")
    else
        println("rcut = $rcut default")
    end

     return rcut
end


fname = "clusterVector.jld2"
if "--fname" in ARGS
    idx = findfirst(x->x=="--fname", ARGS)
    if idx == length(ARGS)
        println("Error: --fname requires an argument")
        exit(1)
    end
    fname = ARGS[idx+1]

end

clusterVector = DC.jldopen(fname)["clusterVector"]
refCNA = clusterVector.vec[1].CNA
if "--ref" in ARGS
    idx = findfirst(x->x=="--ref", ARGS)
    if idx == length(ARGS)
        println("Error: --ref requires an argument")
        exit(1)
    end
    refstructure = DC.read_xyz(ARGS[idx+1])
    rcut = get_rcut()
    refCNA = DC.getCNAProfile(refstructure, rcut)
end

refCNA2 = clusterVector.vec[2].CNA
if "--ref2" in ARGS
    idx = findfirst(x->x=="--ref2", ARGS)
    if idx == length(ARGS)
        println("Error: --ref2 requires an argument")
        exit(1)
    end
    refstructure = DC.read_xyz(ARGS[idx+1])
    rcut = get_rcut()
    refCNA2 = DC.getCNAProfile(refstructure, rcut)
end


sims, energies = DC.getSimsAndEnergies(clusterVector, refCNA)
sims2, energies = DC.getSimsAndEnergies(clusterVector, refCNA2)


toLabel = Vector{Int64}(undef, 0)
if "--label" in ARGS
    idx = findfirst(x->x=="--label", ARGS)
    lstart = 1
    lend = 0
    if occursin(",", ARGS[idx+1])
        lrange = split(ARGS[idx+1], ",")
        lstart = parse(Int64, lrange[1])
        lend = parse(Int64, lrange[2])
    else
        lend = parse(Int64, ARGS[idx+1])
    end
    toLabel = [i for i in lstart:lend]
end

if "-s" in ARGS || "--schebarchov" in ARGS
    u = ["ICO", "ANTIICO", "TWICO", "TRICO", "POLYICO", "EXICO", "DEC", "TWI", "FCC", "AMB"]
    cols = ["red", "purple", "maroon", "maroon", "maroon", "slateblue", "grey", "orange", "blue", "pink"]
    Line2D = PyPlot.matplotlib.lines.Line2D
    legendelements = [Line2D([0], [0], marker="o", label=u[x], color=cols[x], linestyle="") for x in 1:length(u)]
    sc = Vector{String}(undef, clusterVector.N[])
    sidx = findfirst(x->x=="-s" || x=="--schebarchov", ARGS)
    onlyplot = nothing
    if sidx+1 <= length(ARGS) && ARGS[sidx+1] in u
        onlyplot = ARGS[sidx+1]
    end
        

    sclass = readlines(open("schebarchovClass.txt"))
    for i in 1:length(sclass)
        sclass[i] = split(sclass[i])[1]
    end
    for i in 1:length(sc)
        sc[i] = cols[findfirst(x->x==sclass[i], u)]
    end
    
    toplot = nothing
    c = sc
    if onlyplot != nothing
        toplot = findall(x->x==onlyplot, sclass)
        c = "k"
    end
    

    if "--ref2" in ARGS
        scatter(sims, sims2, c=c, s=1)
        if toplot != nothing
            scatter(sims[toplot], sims2[toplot], c=sc[toplot], s=1)
        end
        for i in toLabel
            annotate(string(i), (sims[i], sims2[i]), fontsize=8)
        end
    else

        scatter(sims, energies, c=c, s=1)
        if toplot != nothing
            scatter(sims[toplot], energies[toplot], c=sc[toplot], s=1)
        end
        for i in toLabel
            annotate(string(i), (sims[i], energies[i]), fontsize=8)
        end
    end
    xlim(0,1)
    legend(handles=legendelements, loc="upper right", fontsize=8)
elseif "-g" in ARGS || "--gmm" in ARGS
    c::Vector{Int64} = [parse(Int64, x) for x in readlines(open("gmmclasses.txt"))]
    cmap = "tab10"
    if length(unique(c)) > 10
        cmap = "tab20"
    end
    if "--ref2" in ARGS
        scatter(sims, sims2, c=c, s=1, cmap=cmap)
    else
        scatter(sims, energies, c=c, s=1, cmap=cmap)
    end
    xlim(0,1)
else
    if "--ref2" in ARGS
        scatter(sims, sims2, s=1)
    else
        scatter(sims, energies, s=1)
    end
    xlim(0,1)

end

show()

