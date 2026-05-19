#!/usr/bin/env julia

if "--help" in ARGS || "-h" in ARGS
    println()
    println("Usage: classify_clusterVector.jl [options]")
    println("Options:")
    println("  --fname, -f <filename>   JLD2 file containing the clusterVector (default: clusterVector.jld2)")
    println("  --rcut, -r <value>      Cutoff distance for classification (default: 1.0)")
    println()
    exit(0)
end

using DAC
import .DAC
const DC = DAC

fname = "clusterVector.jld2"
if "--fname" in ARGS || "-f" in ARGS
    idx = findfirst(x -> x == "--fname" || x == "-f", ARGS)
    if idx == length(ARGS)
        println("Error: --fname requires an argument")
        exit(1)
    end
    fname = ARGS[idx+1]

end


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


clusterVector = DC.jldopen(fname)["clusterVector"]
refCNA = clusterVector.vec[1].CNA
classes = Vector{String}(undef, clusterVector.N[])
sims = Vector{Float64}(undef, clusterVector.N[])
n = clusterVector.N[]
ndone = Threads.Atomic{Int64}(0)
println(ndone)
Threads.@threads for i in 1:n
    if ndone[] % 100 == 0
        print("\r $(ndone[]) / $n")
    end
    classes[i] = DC.classifyCluster(clusterVector.vec[i].positions, rcut)
    sims[i] = DC.getCNASimilarity(clusterVector.vec[i].CNA, refCNA)
    Threads.atomic_add!(ndone, 1)
end

schebarchovclassfile = open("schebarchovClass.txt", "w")
for i in 1:length(classes)
    println(schebarchovclassfile, "$(classes[i]) $(clusterVector.vec[i].energy) $(sims[i])")
end
close(schebarchovclassfile)

println()
