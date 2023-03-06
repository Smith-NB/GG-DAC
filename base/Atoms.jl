abstract type Calculator end


"""
    Cluster

# Arguments

- `forrmula::Dict{String, Int64}`: The chemical formula of the cluster, stored as a element keyed Dict type.
- `positions::Array{Float64}`: The Euclidian coordinates of all atoms in the cluster.
- `cell::Array{Float64}`: The dimensions of the lattice cell.
- `calculator::Calculator`: The calculator used to calculate energies, forces, etc.
- `energy::Float64`: Optional. The energy of the cluster. Typically calculated by the calculator.
- `energies::Array{Float64}`: Optional. The energy of each atom in the cluster. Typically calculated by the calculator.
- `forces::Array{Float64}`: Optional. The forces acting on each atom. Typically calculated by the calculator.
- `stresses::Array{Float64}: Optional. The stresses experienced by each atom in the cluster. Typically calculated by the calculator.
"""
mutable struct Cluster
	formula::Dict{String, Int64}
	positions::Array{Float64}
	cell::Array{Float64}
	calculator::Calculator
	energy::Float64
	energies::Array{Float64}
	forces::Array{Float64}
	stresses::Array{Float64}
	Cluster(formula::Dict{String, Int64}, 
			positions::Array{Float64}, 
			cell::Array{Float64}
			) = new(formula, positions, cell)
end

Cluster(formula) = Cluster(formula, zeros(Float64, sum(get.([formula], keys(formula), nothing)), 3), zeros(Float64, 3, 3), nothing)
Cluster(formula, positions) = Cluster(formula, positions, zeros(Float64, 3, 3), nothing)
Cluster(formula, positions, cell) = Cluster(formula, positions, cell, nothing)
	
"""
	getNAtoms(atoms::Cluster)

Returns the atom count of a Cluster type.
"""
function getNAtoms(atoms::Cluster)
	return trunc(Int, length(atoms.positions)/3)
end

"""
	getDistances(Atoms::Cluster)

Returns the distances between all atoms of a Cluster type.
"""
function getDistances(atoms::Cluster)
	natoms = getNAtoms(atoms)

	d2 = zeros(Float64, natoms, natoms, 3)
	r = zeros(Float64, natoms, natoms)
	
	for i in 1:3
		for j in 1:natoms
			d2[j, :, i] = (atoms.positions[j, i] .- atoms.positions[:, i]).^2
		end
	end
	for ii in 1:natoms
		r[:, ii] = sum(d2[ii, :, :].^2, dims=2).^0.5
	end
	
	#=
	for ii in 1:natoms
		for c in 1:3
				distance_vectors[ii, :, c] = atoms.positions[ii, c] .- atoms.positions[:, c]
		end
		r[:, ii] = sum(distance_vectors[ii, :, :].^2, dims=2).^0.5
	end=#
	return r
end


"""
	read_xyz(filename::String)

Takes a path to a ".xyz" file and loads it as a Cluster type,
which is returned.
"""
function read_xyz(filename::String)
	lines = readlines(filename)
	natoms = parse(Int64, lines[1])
	formula = Dict{String, Int64}()
	positions = zeros(Float64, natoms, 3)
	
	cell = split(split(lines[2], "\"", keepempty=false)[2], " ")
	cell = parse.(Float64, cell)

	for i in 3:length(lines)
		data = split(lines[i], " ", keepempty=false)
		if haskey(formula, data[1])
			formula[data[1]] += 1
		else
			formula[data[1]] = 1
		end
		for j in 2:4
			positions[i-2, j-1] = parse(Float64, data[j])
		end
	end
	return Cluster(formula, positions, cell)
end

"""
	write_xyz(filename::String, atoms::Cluster)

Takes a path (should end in ".xyz") to a file to write and saves a Cluster type
to that path as a .xyz file.
"""
function write_xyz(filename::String, atoms::Cluster)
	
	newfile = open(filename, "w")
	natoms = getNAtoms(atoms)

	print(newfile, natoms, "\n")

	print(newfile, "Lattice=\"", atoms.cell[1])
	for i in 2:9
		print(newfile, " ", atoms.cell[i])
	end
	print(newfile, "\" Properties=species:S:1:pos:R:3 pbc=\"F F F\"\n")

	elements = keys(atoms.formula)
	tab = "        " #8 spaces
	for (element, freq) in atoms.formula
		for i in 1:freq
			print(newfile, element) 
			for j in 1:3
				s = rpad(round(atoms.positions[i, j], sigdigits=10), 11, "0") #pad
				print(newfile, tab, s)
			end
			print(newfile, "\n")
		end
	end
	close(newfile)
end




