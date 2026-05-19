#!/usr/bin/env julia

if "--help" in ARGS || length(ARGS) < 2
    println()
    println("Usage: combine_clusterVectors.jl <start> <final> [--energycutoff <value>] [--fout <filename>]")
    println("Example: combine_clusterVectors.jl 1 10")
    println("Combines cluster vectors from Trial<start> to Trial<final> into a single cluster vector, called totalClusterVector.jld2")
    println("Optional arguments:")
    println("  --energycutoff <value>  Only include clusters with energy less than the specified value")
    println("  --fout <filename>      Specify the output filename (default: totalClusterVector.jld2)")
    println()
    exit(0)
end

using DAC
import .DAC
const DC = DAC

start = parse(Int, ARGS[1])
final = parse(Int, ARGS[2])

clusterVector = DC.jldopen("Trial$start/clusterVector.jld2")["clusterVector"]

dp = 2
ecut = -Inf
if "--energycutoff" in ARGS
    idx = findfirst(==( "--energycutoff"), ARGS)
    ecut = parse(Float64, ARGS[idx + 1])
end
for i in 1:clusterVector.N[]
    if clusterVector.vec[i].energy > ecut
        break
    end
    e = round(clusterVector.vec[i].energy, digits=dp)
    clusterVector.vec[i] = DC.ClusterCompressed(clusterVector.vec[i].positions, e, clusterVector.vec[i].CNA, clusterVector.vec[i].ID)
end

for i in start:final
toAdd = DC.jldopen("Trial$i/clusterVector.jld2")["clusterVector"]
    for j in 1:toAdd.N[]
        DC.addToVector!(toAdd.vec[j], clusterVector, dp)
    end
    print("\radded $i")
end
println()

fout = "totalClusterVector.jld2"
if "--fout" in ARGS
    idx = findfirst(==( "--fout"), ARGS)
    fout = ARGS[idx + 1]
else
    println("No output file specified, using default: totalClusterVector.jld2")
end

DC.@save fout clusterVector
