#!/usr/bin/env julia

if "--help" in ARGS
    println("Usage: getxyzofLESs.jl [clusterVector.jld2] [num]")
    println("  clusterVector.jld2: The JLD2 file containing the cluster vector (default: clusterVector.jld2)")
    println("  num: The number of LESs to process (default: 100)")
    exit(0)
end


using DAC
import .DAC
const DC = DAC


fname = "clusterVector.jld2"
num = 100
if length(ARGS) > 0
    fname = ARGS[1]
end
if length(ARGS) > 1
    num = parse(Int, ARGS[2])
end

if !isdir("LES")
    mkdir("LES")
end

clusterVector = DC.jldopen(fname)["clusterVector"]
l = length("$num")
formula = Dict("Cu" => 78)
for i in 1:num
    pos = clusterVector.vec[i].positions
    xyzfname = "$i"
    while length(xyzfname) < l
        xyzfname = "0$xyzfname"
    end
    xyzfname = "$xyzfname.xyz"
    DC.write_xyz("LES/$xyzfname", pos, formula, 15.0)
end
