abstract type Atoms end
abstract type Calculator end
abstract type _py_Calculator <: Calculator end


struct ClusterCompressed
	positions::Array{Float64}
	energy::Float64
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	ID::Int64
end

mutable struct ClusterVector
	vec::Vector{ClusterCompressed}
	N::Threads.Atomic{Int64}
	lock::ReentrantLock
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
	distances::Matrix{Float64}
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
	distances::Matrix{Float64},
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}},
	validCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool
			) = new(formula, positions, cell, energy, energies, forces, stresses, distances, CNA, validCNA, validEnergies, validForces, validStresses)
	Cluster(formula::Dict{String, Int64},
	positions::Matrix{Float64},
	cell::Matrix{Float64},
	energy::Float64,
	energies::Vector{Float64},
	forces::Matrix{Float64},
	stresses::Array{Float64},
	distances::Matrix{Float64},
	CNA::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}},
	validCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool,
	calculator::Calculator
			) = new(formula, positions, cell, energy, energies, forces, stresses, distances, CNA, validCNA, validEnergies, validForces, validStresses, calculator)

end

#Cluster(formula::Dict{String, Int64}) = Cluster(formula, zeros(Float64, sum(get.([formula], keys(formula), nothing)), 3), zeros(Float64, 3, 3))
function Cluster(formula::Dict{String, Int64}) 
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(
		formula, zeros(Float64, N, 3), zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), Matrix{Float64}(undef, natoms, natoms),
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(), 
		false, false, false, false
		)
end 

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(formula, positions, zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), getDistances(positions),
		Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(),
		false, false, false, false
		)
end

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64}, cell::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))

	Cluster(formula, positions, cell, 
		0.0, zeros(Float64, N), zeros(Float64, 3, N), 
		zeros(Float64, 3, 3, N), getDistances(positions),
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
function getCNAProfile(atoms::Cluster) 
	if atoms.validCNA
		return atoms.CNA
	else
		return nothing
	end
end

setFormula(atoms::Cluster, formula::Dict{String, Int64}) = atoms.formula = formula
function setPositions!(atoms::Cluster, positions::Matrix{Float64})
	atoms.positions = positions
	atoms.validCNA = false
	atoms.validEnergies = false
	atoms.validForces = false
	atoms.validStresses = false
end
function moveAtoms!(atoms::Cluster, dr::LinearAlgebra.Adjoint{Float64, Matrix{Float64}})
	atoms.positions += dr
end
function moveAtoms!(atoms::Cluster, dr::Matrix{Float64})
	atoms.positions += dr
end
setCell!(atoms::Cluster, cell::Matrix{Float64}) = atoms.cell = cell
setCalculator!(atoms::Cluster, calculator::Calculator) = atoms.calculator = calculator
function setEnergies!(atoms::Cluster, energies::Vector{Float64}) 
	atoms.energies = energies
	atoms.energy = sum(energies)
	validEnergies = true
end
setEnergy!(atoms::Cluster, energy::Float64) = atoms.energy = energy
setValidCNA!(atoms::Cluster, valid::Bool) = atoms.validCNA = valid
setValidEnergies!(atoms::Cluster, valid::Bool) = atoms.validEnergies = valid
setValidForces!(atoms::Cluster, valid::Bool) = atoms.validForces = valid
setValidStresses!(atoms::Cluster, valid::Bool) = atoms.validStresses = valid
setCalculator!(atoms::Cluster, calc::Calculator) = atoms.calculator = calc
function setCNAProfile!(atoms::Cluster, rcut::Float64)
	atoms.CNA = getCNAProfile(atoms, rcut)
	atoms.validCNA = true
end
setCNAProfile!(atoms::Cluster, cna::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}) = atoms.CNA, atoms.validCNA = cna, true


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
    natoms = getNAtoms(coordinates)


   r = Matrix{Float64}(undef, natoms, natoms)
    for i in 1:natoms
        r[i, i] = 0
        for j in i+1:natoms
        	d2_1 = (coordinates[i, 1] - coordinates[j, 1])^2
        	d2_2 = (coordinates[i, 2] - coordinates[j, 2])^2
        	d2_3 = (coordinates[i, 3] - coordinates[j, 3])^2

            r[j, i] = (d2_1 + d2_2 + d2_3)^0.5
            r[i, j] = r[j, i]
        end
    end

    return r
end

"""
	getDistances(Atoms::Cluster)

Returns the distances between all atoms of a Cluster type.
"""
function getDistances(atoms::Cluster)
	return atoms.distances
	#return getDistances(atoms.positions)
end

"""
	getDistances!(Atoms::Cluster)

sets the atoms.distances attribute to the distances between all atoms.
"""
function setDistances!(atoms::Cluster)
    natoms = getNAtoms(atoms)

    for i in 1:natoms
        atoms.distances[i, i] = 0
        for j in i+1:natoms
        	d2_1 = (atoms.positions[i, 1] - atoms.positions[j, 1])^2
        	d2_2 = (atoms.positions[i, 2] - atoms.positions[j, 2])^2
        	d2_3 = (atoms.positions[i, 3] - atoms.positions[j, 3])^2

            atoms.distances[j, i] = (d2_1 + d2_2 + d2_3)^0.5
            atoms.distances[i, j] = atoms.distances[j, i]
        end
    end

    return nothing
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
	inclusionRadiusOfCluster(coords::Matrix{Float64})

Returns the radous of a sphere that can completely inclose the cluster with a radium from the center of mass of
	to the most out atom from the centre of mass.
"""
function getInclusionRadiusOfCluster(coords::Matrix{Float64})
    maxSize = -Inf
    N = getNAtoms(coords)
    for i in 1:N
        for j in i+1:N
            xDist = coords[i, 1] - coords[j, 1]
            yDist = coords[i, 2] - coords[j, 2]
            zDist = coords[i, 3] - coords[j, 3]
            dist = sqrt(xDist^2 + yDist^2 + zDist^2)

            if dist > maxSize
                maxSize = dist
            end
        end
    end
    #=
    if maxSize == -Inf || maxSize == Inf
    	print("INFALERT______________")
		for i in 1:N
	        for j in i+1:N
	            xDist = coords[i, 1] - coords[j, 1]
	            yDist = coords[i, 2] - coords[j, 2]
	            zDist = coords[i, 3] - coords[j, 3]
	            dist = sqrt(xDist^2 + yDist^2 + zDist^2)
	            println("MAXSIZE $i $j $dist $(coords[i, :]) $(coords[j, :]) $xDist $yDist $zDist")
	            if dist > maxSize
	                maxSize = dist
	            end
	        end
	    end
    end
	=#
    return maxSize/2.0
end	

function getInclusionRadiusOfCluster(atoms::Cluster)
    return getInclusionRadiusOfCluster(atoms.positions)
end


"""
	getCentreOfCluster(coordinates::Matrix{Float64})

Returns the center coordinate of a coordinate Matrix by calculating the center of mass,
	assuming all atoms have equal mass.
"""
function getCentreOfCluster(coordinates::Matrix{Float64})
	N = getNAtoms(coordinates)
	centre = zeros(1, 3)
	for i in 1:3
		for j in 1:N
			centre[i] += coordinates[j, i]
		end
		centre[i] /= N
	end

	return centre
end

"""
	getCentreOfCluster(coordinates::Matrix{Float64})

Returns the center coordinate of a coordinate Matrix by calculating the center of mass,
	with `masses` given in a 1xN matrix.
"""
function getCentreOfCluster(coordinates::Matrix{Float64}, masses::Matrix{Float64})
	N = getNAtoms(coordinates)
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

function classifyAtoms(coordinates::Matrix{Float64}, rcut::Float64)

	nCNA = getNormalCNAProfile(coordinates, rcut)
	N = getNAtoms(coordinates)
	atomClass = Vector{Int64}(undef, N)
	sigs = [(4, 2, 1), (3, 1, 1), (3, 2, 2), (2, 1, 1), (4, 2, 2), (5, 5, 5), (2, 0, 0), (4, 3, 3)]

	for n in 1:N
		# Set non-present sigs to 0.
		for s in sigs
			if !haskey(nCNA[n], s)
				nCNA[n][s] = 0
			end
		end

		# Classify atom.
		if nCNA[n][(5, 5, 5)] >= 6
			atomClass[n] = 11
		elseif nCNA[n][(5, 5, 5)] >= 3
			atomClass[n] = 13
		elseif nCNA[n][(5, 5, 5)] == 2
			atomClass[n] = 6
		elseif nCNA[n][(5, 5, 5)] == 1
			if nCNA[n][(2, 0, 0)] == 2
				atomClass[n] = 12
			else
				atomClass[n] = 7
			end
		elseif nCNA[n][(4, 2, 2)] >= 5
			atomClass[n] = 5
		elseif nCNA[n][(4, 2, 1)] >= 4
			atomClass[n] = 1
		elseif nCNA[n][(3, 1, 1)] >= 3 && nCNA[n][(3, 2, 2)] == 0
			atomClass[n] = 2
		elseif nCNA[n][(2, 1, 1)] >= 3
			atomClass[n] = 3
		elseif nCNA[n][(2, 1, 1)] == 2 && nCNA[n][(3, 1, 1)] == 2
			atomClass[n] = 4
		elseif nCNA[n][(3, 1, 1)] >= 3 && nCNA[n][(3, 2, 2)] >= 1
			atomClass[n] = 8
		elseif nCNA[n][(3, 1, 1)] >= 1 && nCNA[n][(2, 0, 0)] == 2 && nCNA[n][(2, 1, 1)] >= 1
			atomClass[n] = 9
		elseif nCNA[n][(4, 3, 3)] >= 1 && nCNA[n][(2, 0, 0)] == 2
			atomClass[n] = 14
		else
			atomClass[n] = 10
		end
	end

	return atomClass
end

classifyAtoms(cluster::Cluster, rcut::Float64) = classifyAtoms(cluster.positions, rcut)

function classifyCluster(coordinates::Matrix{Float64}, rcut::Float64)

	tag = classifyAtoms(coordinates, rcut)

	icoAtoms = Vector{Int64}()
	icoBonds = Vector{Int64}()
	ico12Atoms = Vector{Int64}()
	icoCores = Vector{Int64}()
	nIcoAtoms = 0
	nIcoBonds = 0
	nIco12Atoms = 0
	nIco12Bonds = 0
	nIcoCores = 0

	nFCC = 0
	nHCP = 0
	nIcoSpines = 0
	nPartIcoCores = 0
	nAmb = 0
	nClass14 = 0

	N = getNAtoms(coordinates)


	# Count atoms of each type and store some atom indicies for Ih cases.
	for n in 1:N
		if tag[n] in [6, 7, 11, 12, 13]
			push!(icoAtoms, n)
			nIcoAtoms += 1
			push!(icoBonds, 0)
			if tag[n] == 11
				nIcoCores += 1
				push!(icoCores, n)
			elseif tag[n] == 12
				push!(ico12Atoms, n)
				nIco12Atoms += 1
			elseif tag[n] == 13
				nPartIcoCores += 1
			end
		elseif tag[n] in [1, 2, 3, 4]
			nFCC += 1
		elseif tag[n] == 5
			nHCP += 1
		elseif tag[n] == 10
			nAmb += 1
		elseif tag[n] == 14
			nClass14 += 1
			push!(ico12Atoms, n)
			nIco12Atoms += 1
		end

	end

	# Determine which ico spine atoms are bonding.
	for i in 1:nIcoAtoms
		for j in i+1:nIcoAtoms
			#calculate distance between soube atoms
			d = ((coordinates[icoAtoms[i], 1] - coordinates[icoAtoms[j], 1])^2
				+(coordinates[icoAtoms[i], 2] - coordinates[icoAtoms[j], 2])^2
				+(coordinates[icoAtoms[i], 3] - coordinates[icoAtoms[j], 3])^2
				)^0.5
			if d <= rcut
				icoBonds[i] += 1
				icoBonds[j] += 1
				nIcoBonds += 1
			end
		end		
	end

	# Check which 12/14 atoms are bonding (antiMackay 555/433)
	for i in 1:nIco12Atoms
		for j in i+1:nIco12Atoms
			d = ((coordinates[ico12Atoms[i], 1] - coordinates[ico12Atoms[j], 1])^2
				+(coordinates[ico12Atoms[i], 2] - coordinates[ico12Atoms[j], 2])^2
				+(coordinates[ico12Atoms[i], 3] - coordinates[ico12Atoms[j], 3])^2
				)^0.5
			if d <= rcut
				nIco12Bonds += 1
			end
		end	
	end

	# Check if all cores are true cores. (a core bonding to a 12 or 14 atom is a false core.)
	for i in 1:nIcoCores
		for j in 1:nIco12Atoms
			d = ((coordinates[icoCores[i], 1] - coordinates[ico12Atoms[j], 1])^2
				+(coordinates[icoCores[i], 2] - coordinates[ico12Atoms[j], 2])^2
				+(coordinates[icoCores[i], 3] - coordinates[ico12Atoms[j], 3])^2
				)^0.5
			if d <= rcut
				if nIcoCores > 1
					nIcoCores -= 1
				else
					nIcoCores = -1
				end
				break
			end
		end
	end

	# Determine and return class.
	if nIcoCores != 0
		if nIcoCores > 3
			return "POLYICO"
		elseif nIcoCores > 2
			return "TRICO"
		elseif nIcoCores > 1
			return "TWICO"
		elseif nIcoCores == -1
			return "EXICO"
		else
			return nIco12Bonds != 0 ? "ANTIICO" : "ICO"
		end
	elseif nIcoAtoms != 0 && nPartIcoCores == 0
		return "DEC"
	elseif nIcoAtoms == 0 && nHCP + nFCC > nAmb
		if nHCP == 0
			return "FCC"
		elseif nFCC == 0
			return "HCP"
		end
		return "TWI"
	else
		return "AMB"
	end
end

classifyCluster(cluster::Cluster, rcut::Float64) = classifyCluster(cluster.positions, rcut)

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

"""
	write_xyz(filename::String, positions::Matrix{Float64})

Takes a path (should end in ".xyz") to a file to write and saves a Cluster type
to that path as a .xyz file.
"""
function write_xyz(filename::String, positions::Matrix{Float64}, formula::Dict{String, Int64}, cell::Float64)
	
	newfile = open(filename, "w")
	natoms = getNAtoms(positions)

	print(newfile, natoms, "\n")

	print(newfile, "Lattice=\"$(cell) 0.0 0.0 0.0 $(cell) 0.0 0.0 0.0 $(cell)")
	print(newfile, "\" Properties=species:S:1:pos:R:3 pbc=\"F F F\"\n")

	elements = keys(formula)
	tab = "        " #8 spaces
	for (element, freq) in formula
		for i in 1:freq
			print(newfile, element) 
			for j in 1:3
				s = rpad(round(positions[i, j], sigdigits=10), 11, "0") #pad
				print(newfile, tab, s)
			end
			print(newfile, "\n")
		end
	end
	close(newfile)
end

"""
	write_xyz(filename::String, atoms::Cluster)

Takes a path (should end in ".xyz") to a file to write and saves a Cluster type
to that path as a .xyz file. Also stores the `tags' of each atom.
"""
function write_xyz(filename::String, atoms::Cluster, tags::Vector{Int64})
	newfile = open(filename, "w")
	natoms = getNAtoms(atoms)

	print(newfile, natoms, "\n")

	print(newfile, "Lattice=\"", atoms.cell[1])
	for i in 2:9
		print(newfile, " ", atoms.cell[i])
	end
	print(newfile, "\" Properties=species:S:1:pos:R:3:tags:I:1 pbc=\"F F F\"\n")

	elements = keys(atoms.formula)
	tab = "        " #8 spaces
	for (element, freq) in atoms.formula
		for i in 1:freq
			print(newfile, element) 
			for j in 1:3
				s = rpad(round(atoms.positions[i, j], sigdigits=10), 11, "0") #pad
				print(newfile, tab, s)
			end
			print(newfile, tab * string(tags[i]) * "\n")
		end
	end
	close(newfile)
end

function formulaDictToString(formula::Dict{String, Int64})
	s = ""
	for key in keys(formula)
		s *= string(key)
		s *= string(formula[key])
	end

	return s
end

function formulaDictToString(atoms::Cluster) 
	return formulaDictToString(atoms.formula)
end

"""
	view(atoms::Cluster)

Takes a Cluster type and calls ASE gui module for visualisation.
"""
function view(atoms::Cluster)
	_view = pyimport("ase.visualize").view
	_atoms = pyimport("ase").Atoms

	atoms = _atoms(formulaDictToString(atoms.formula), positions = atoms.positions, cell = atoms.cell)
	_view(atoms)
end

function view(coords::Matrix{Float64}, formula::Dict{String, Int64})
	_view = pyimport("ase.visualize").view
	_atoms = pyimport("ase").Atoms

	atoms = _atoms(formulaDictToString(formula), positions = coords)
	_view(atoms)
end

function view(atoms::ClusterCompressed, formula::Dict{String, Int64})
	view(atoms.positions, formula)

end

