#!/usr/bin/env julia

if "--help" in ARGS
    println("Usage: tagxyzs.jl [directory] [rcut]")
    println("Tags atoms in .xyz files in the specified directory using Schebarchov's method.")
    println("  directory: directory containing .xyz files (default: LES)")
    println("  rcut: cutoff distance for classifying atoms (default: 1.0 or from BHA_run_script.jl or DACA_run_file.jl if found)")
    exit(0)
end

using DAC
import .DAC
const DC = DAC



dname = "LES"
if length(ARGS) != 0
    dname = ARGS[1]
end

rcut = 1.0
runfilename = nothing
if isfile("BHA_run_script.jl")
    runfilename = "BHA_run_script.jl"
elseif isfile("DACA_run_file.jl")
    runfilename = "DACA_run_file.jl"
end

if length(ARGS) > 1
    rcut = parse(Float64, ARGS[2])
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


for f in readdir(dname)
    if !endswith(f, ".xyz") continue end
    atoms = DC.read_xyz("$dname/$f")
    tags = DC.classifyAtomsSchebarchov(atoms, rcut)
    DC.write_xyz("$dname/$f", atoms, tags)
end

