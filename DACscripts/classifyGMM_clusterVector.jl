#!/usr/bin/env julia

if "--help" in ARGS || "-h" in ARGS
    println()
    println("Usage: classifyGMM_clusterVector.jl [options]")
    println("Options:")
    println("  --fname, -f <filename>   Specify the JLD2 file containing the cluster vector (default: clusterVector.jld2)")
    println("  --rcut, -r <value>      Specify the cutoff radius for classifying structures (default: 1.0)")
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

gmmfname = "GMM_PCA.jld2"
gmm = DC.jldopen(gmmfname)["gmm"]

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
classes = Vector{String}(undef, clusterVector.N[])
n = clusterVector.N[]
ndone = Threads.Atomic{Int64}(0)
println(ndone)
classMatrix::Matrix{Float64} = DC.getClassMatrix(clusterVector, rcut, DC.getClasses())
gmmclasses = DC.getStructureClasses(clusterVector, gmm, classMatrix'[:,:])
gmmclassfile = open("gmmclasses.txt", "w")
for i in 1:length(gmmclasses)
    println(gmmclassfile, gmmclasses[i])
end
