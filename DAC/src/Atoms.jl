abstract type Atoms end
abstract type Calculator end
abstract type _py_Calculator <: Calculator end


struct ClusterCompressed
	positions::Array{Float64}
	energy::Float64
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
end

mutable struct ClusterVector
	vec::Vector{ClusterCompressed}
	N::Int64
	lock::Bool
end

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
mutable struct Cluster <: Atoms
	formula::Dict{String, Int64}
	positions::Matrix{Float64}
	cell::Matrix{Float64}
	energy::Float64
	energies::Vector{Float64}
	forces::Matrix{Float64}
	stresses::Array{Float64}
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	validCNA::Bool
	validEnergies::Bool
	validForces::Bool
	validStresses::Bool
	calculator::Calculator
	Cluster(formula::Dict{String, Int64},
	positions::Matrix{Float64},
	cell::Matrix{Float64},
	energy::Float64,
	energies::Vector{Float64},
	forces::Matrix{Float64},
	stresses::Array{Float64},
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}},
	validCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool
			) = new(formula, positions, cell, energy, energies, forces, stresses, CNA, validCNA, validEnergies, validForces, validStresses)
	Cluster(formula::Dict{String, Int64},
	positions::Matrix{Float64},
	cell::Matrix{Float64},
	energy::Float64,
	energies::Vector{Float64},
	forces::Matrix{Float64},
	stresses::Array{Float64},
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}},
	validCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool,
	calculator::Calculator
			) = new(formula, positions, cell, energy, energies, forces, stresses, CNA, validCNA, validEnergies, validForces, validStresses, calculator)

end

#Cluster(formula::Dict{String, Int64}) = Cluster(formula, zeros(Float64, sum(get.([formula], keys(formula), nothing)), 3), zeros(Float64, 3, 3))
function Cluster(formula::Dict{String, Int64}) 
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(
		formula, zeros(Float64, N, 3), zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), 
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(), 
		false, false, false, false
		)
end 

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(formula, positions, zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), 
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(),
		false, false, false, false
		)
end

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64}, cell::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(formula, positions, cell, 
		0.0, zeros(Float64, N), zeros(Float64, 3, N), 
		zeros(Float64, 3, 3, N), 
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(),
		false, false, false, false
		)
	
	#=
	Cluster(formula, positions, cell, 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), 
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(),
		false, false
		)
	=#
end


Base.copy(c::Cluster) = Cluster(c.formula, c.positions, c.cell, c.energy, c.energies, c.forces, c.stresses, 
							c.CNA, c.validCNA, c.validEnergies, c.validForces, c.validStresses, c.calculator)

getFormula(atoms::Cluster) = atoms.formula
getPositions(atoms::Cluster) = atoms.positions
getCell(atoms::Cluster) = atoms.cell
getEnergy(atoms::Cluster) = atoms.energy
getEnergies(atoms::Cluster) = atoms.energies
function getForces!(atoms::Cluster) 
	calculateForces!(atoms, getCalculator(atoms))
	return atoms.forces
end
getForces(atoms::Cluster) = atoms.forces
getStresses(atoms::Cluster) = atoms.stresses
getCNA(atoms::Cluster) = atoms.CNA
getValidCNA(atoms::Cluster) = atoms.validCNA
getValidEnergies(atoms::Cluster) = atoms.validEnergies
getValidForces(atoms::Cluster) = atoms.validForces
getValidStresses(atoms::Cluster) = atoms.validStresses
getCalculator(atoms::Cluster) = atoms.calculator


setFormula(atoms::Cluster, formula::Dict{String, Int64}) = atoms.formula = formula
function setPositions!(atoms::Cluster, positions::Matrix{Float64})
	atoms.positions = positions
	atoms.validCNA = false
	atoms.validEnergies = false
	atoms.validForces = false
	atoms.validStresses = false
end
setCell!(atoms::Cluster, cell::Matrix{Float64}) = atoms.cell = cell
setCalculator!(atoms::Cluster, calculator::Calculator) = atoms.calculator = calculator
function setEnergies!(atoms::Cluster, energies::Vector{Float64}) 
	atoms.energies = energies
	atoms.energy = sum(energies)
	validEnergies = true
end
setValidCNA!(atoms::Cluster, valid::Bool) = atoms.validCNA = valid
setValidEnergies!(atoms::Cluster, valid::Bool) = atoms.validEnergies = valid
setValidForces!(atoms::Cluster, valid::Bool) = atoms.validForces = valid
setValidStresses!(atoms::Cluster, valid::Bool) = atoms.validStresses = valid
setCalculator!(atoms::Cluster, calc::Calculator) = atoms.calculator = calculator
"""
	hasCNAProfile(atoms::Cluster)

Returns boolean for if the CNA field of a Cluster instance has been set.
"""
function hasCNAProfile(atoms::Cluster)
	try
		if atoms.CNA != nothing
			return true
		else
			return false
		end
	catch UndefRefError
		return false
	end
end

"""
	getNAtoms(coordinates::Matrix{Float64})

Returns the atom count of a Matrix type.
"""
getNAtoms(coordinates::Matrix{Float64}) = trunc(Int, length(coordinates)/3)

"""
	getNAtoms(atoms::Cluster)

Returns the atom count of a Cluster type.
"""
getNAtoms(atoms::Cluster) = getNAtoms(atoms.positions)


"""
	getDistances(coordinates::Matrix{Float64})

Returns the distances between all atoms from a Maxtrix of atomic coordinates.
"""
function getDistances(coordinates::Matrix{Float64})
	natoms = trunc(Int, length(coordinates)/3)


	d2 = zeros(Float64, natoms, natoms, 3)
	r = zeros(Float64, natoms, natoms)
	
	for i in 1:3
		for j in 1:natoms
			d2[j, :, i] = (coordinates[j, i] .- coordinates[:, i]).^2
		end
	end
	for ii in 1:natoms
		r[:, ii] = sum(d2[ii, :, :], dims=2).^0.5
	end
	
	return r
end

"""
	getDistances(Atoms::Cluster)

Returns the distances between all atoms of a Cluster type.
"""
function getDistances(atoms::Cluster)
	return getDistances(atoms.positions)
end

"""
	getNeighboursList(coordinates::Matrix{Float64}, maxBondingDistance::Number)

Returns the neighbours for all atoms from a coordinate Matrix, neighbours being other atoms
	closer than the maxBongingDistance
"""
function getNeighboursList(coordinates::Matrix{Float64}, maxBondingDistance::Number)
	N = getNAtoms(coordinates)
	neighbourList = [Int64[] for x in 1:N]
	r = getDistances(coordinates)

	for i in 1:N-1
		for j in i+1:N
			if r[j, i] < maxBondingDistance
				push!(neighbourList[i], j)
				push!(neighbourList[j], i)
			end
		end
	end

	return neighbourList

end

"""
	getNeighboursList(atoms::Cluster, maxBondingDistance::Float64)

Returns the neighbours for all atoms from a Cluster type, neighbours being other atoms
	closer than the maxBongingDistance
"""
function getNeighboursList(atoms::Cluster, maxBondingDistance::Float64)
	return getNeighboursList(atoms.positions, maxBondingDistance)
end


"""
	getCentreOfCluster(coordinates::Matrix{Float64})

Returns the center coordinate of a coordinate Matrix by calculating the center of mass,
	assuming all atoms have equal mass.
"""
function getCentreOfCluster(coordinates::Matrix{Float64})
	N = getNAtoms(coordinates)
	masses = fill(1.0, 1, N)
	centre = masses * coordinates ./ N
	return centre
end


"""
	getCentreOfCluster(atoms::Cluster)

Returns the center coordinate of a cluster by calculating the center of mass,
	assuming all atoms have equal mass.
"""
function getCentreOfCluster(atoms::Cluster)
	return getCentreOfCluster(atoms.positions)
end


"""
	centreCluster!(atoms::Cluster)

Centres a cluster within its cell.
"""
function centreCluster!(atoms::Cluster)
	centre = getCentreOfCluster(atoms)
	cellcentre = [atoms.cell[i]/2 for i in [1, 5, 9]]
	translation = cellcentre' - centre
	atoms.positions .+= translation
end

function read_xyzs(filename::String, formula::Dict{String, Int64})
	lines = readlines(filename)
	natoms = parse(Int64, lines[1])
	linesPerCluster = natoms + 2
	nclusters = trunc(Int, length(lines)/(natoms+2))
	clusters = Array{Cluster}(undef, nclusters)

	for i in 1:nclusters
		l = (i-1)*linesPerCluster
		positions = zeros(Float64, natoms, 3)
		cell = reshape(split(split(lines[l+2], "\"", keepempty=false)[2], " "), 3, 3)
		cell = parse.(Float64, cell)
		for j in l+3:l+linesPerCluster
			data = split(lines[j], " ", keepempty=false)
			for k in 2:4
				positions[j-l-2, k-1] = parse(Float64, data[k])
			end
		end
		clusters[i] = Cluster(formula, positions, cell)
	end
	return clusters
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
	
	cell = reshape(split(split(lines[2], "\"", keepempty=false)[2], " "), 3, 3)
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




