abstract type Atoms end
abstract type Calculator end
abstract type _py_Calculator <: Calculator end

const CNAProfile = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
const normalCNAProfile = Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}
const CNASig = Tuple{UInt8, UInt8, UInt8}

struct ClusterCompressed
	positions::Array{Float64}
	energy::Float64
	CNA::CNAProfile
	ID::Int64
end


mutable struct ClusterVector
	vec::Vector{ClusterCompressed}
	N::Threads.Atomic{Int64}
	lock::ReentrantLock
end

mutable struct ClusterVectorWithML
	vec::Vector{ClusterCompressed}
	highEVec::Vector{ClusterCompressed}
	eLim::Float64
	MLData::Matrix{UInt8}
	nMLData::Int64
	idsOfMLLabels::Vector{Vector{Int32}}
	idToIndex::Vector{Int32}
	N::Threads.Atomic{Int64}
	lock::ReentrantLock
end

# old, non-thread-safe version for backwards compatability
#=
mutable struct ClusterVector
	vec::Vector{ClusterCompressed}
	N::Int64
	lock::Bool
end
=#

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
	CNA::CNAProfile
	nCNA::normalCNAProfile
	atomClassCount::Vector{UInt8}
	mlLabel::Int64
	validCNA::Bool
	validnCNA::Bool
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
	CNA::CNAProfile,
	nCNA::normalCNAProfile,
	atomClassCount::Vector{UInt8},
	validCNA::Bool,
	validnCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool
			) = new(formula, positions, cell, energy, energies, forces, stresses, distances, CNA, nCNA, atomClassCount, validCNA, validEnergies, validForces, validStresses)
	Cluster(formula::Dict{String, Int64},
	positions::Matrix{Float64},
	cell::Matrix{Float64},
	energy::Float64,
	energies::Vector{Float64},
	forces::Matrix{Float64},
	stresses::Array{Float64},
	distances::Matrix{Float64},
	CNA::CNAProfile,
	nCNA::normalCNAProfile,
	atomClassCount::Vector{UInt8},
	validCNA::Bool,
	validnCNA::Bool,
	validEnergies::Bool,
	validForces::Bool,
	validStresses::Bool,
	calculator::Calculator
			) = new(formula, positions, cell, energy, energies, forces, stresses, distances, CNA, nCNA, atomClassCount, validCNA, validEnergies, validForces, validStresses, calculator)

end

#Cluster(formula::Dict{String, Int64}) = Cluster(formula, zeros(Float64, sum(get.([formula], keys(formula), nothing)), 3), zeros(Float64, 3, 3))
function Cluster(formula::Dict{String, Int64}) 
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(
		formula, zeros(Float64, N, 3), zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), Matrix{Float64}(undef, natoms, natoms),
		CNAProfile(), normalCNAProfile(), Vector{UInt8}(),
		false, false, false, false, false
		)
end 

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))
	Cluster(formula, positions, zeros(Float64, 3, 3), 
		0.0, zeros(Float64, N), zeros(Float64, N, 3), zeros(Float64, N, 3, 3), getDistances(positions),
		CNAProfile(), normalCNAProfile(), Vector{UInt8}(),
		false, false, false, false, false
		)
end

function Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64}, cell::Matrix{Float64})
	N = sum(get.([formula], keys(formula), nothing))

	Cluster(formula, positions, cell, 
		0.0, zeros(Float64, N), zeros(Float64, 3, N), 
		zeros(Float64, 3, 3, N), getDistances(positions),
		CNAProfile(), normalCNAProfile(), Vector{UInt8}(),
		false, false, false, false, false
		)
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
getNormalCNA(atoms::Cluster) = atoms.nCNA
getValidCNA(atoms::Cluster) = atoms.validCNA
getValidnCNA(atoms::Cluster) = atoms.validnCNA
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
	return nothing
end
function moveAtoms!(atoms::Cluster, dr::LinearAlgebra.Adjoint{Float64, Matrix{Float64}})
	atoms.positions += dr
	return nothing
end
function moveAtoms!(atoms::Cluster, dr::Matrix{Float64})
	atoms.positions += dr
	return nothing
end
setCell!(atoms::Cluster, cell::Matrix{Float64}) = atoms.cell = cell
setCalculator!(atoms::Cluster, calculator::Calculator) = atoms.calculator = calculator
function setEnergies!(atoms::Cluster, energies::Vector{Float64}) 
	atoms.energies = energies
	atoms.energy = sum(energies)
	validEnergies = true
	return nothing
end
setEnergy!(atoms::Cluster, energy::Float64) = atoms.energy = energy
setValidCNA!(atoms::Cluster, valid::Bool) = atoms.validCNA = valid
setValidnCNA!(atoms::Cluster, valid::Bool) = atoms.validnCNA = valid
setValidEnergies!(atoms::Cluster, valid::Bool) = atoms.validEnergies = valid
setValidForces!(atoms::Cluster, valid::Bool) = atoms.validForces = valid
setValidStresses!(atoms::Cluster, valid::Bool) = atoms.validStresses = valid
setCalculator!(atoms::Cluster, calc::Calculator) = atoms.calculator = calc
function setCNAProfile!(atoms::Cluster, rcut::Float64)
	atoms.CNA = getCNAProfile(atoms.positions, rcut)
	atoms.validCNA = true
	return nothing
end
setTotalCNAProfile!(atoms::Cluster, rcut::Float64) = setCNAProfile!(atoms, rcut)
function setNormalCNAProfile!(atoms::Cluster, rcut::Float64)
	atoms.nCNA = getNormalCNAProfile(atoms.positions, rcut)
	atoms.validnCNA = true
end
function setCNAProfiles!(atoms::Cluster, rcut::Float64)
	atoms.CNA, atoms.nCNA = getTotalAndNormalCNAProfile(atoms.positions, rcut)
	atoms.validCNA = true
	atoms.validnCNA = true
	return nothing
end

setCNAProfile!(atoms::Cluster, cna::CNAProfile) = atoms.CNA, atoms.validCNA = cna, true
setCNAProfiles!(atoms::Cluster, cna::CNAProfile, ncna::normalCNAProfile) = atoms.CNA, atoms.nCNA, atoms.validCNA = cna, ncna, true


setAtomClassCount(atoms::Cluster, classes::normalCNAProfile, nClasses::Int64, returnType::DataType) = atoms.atomClassCount = getFrequencyClassVector(getAtomClasses(atoms.nCNA, classes), nClasses, returnType)

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

function classifyAtomsSchebarchov(coordinates::Matrix{Float64}, rcut::Float64)

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

classifyAtomsSchebarchov(cluster::Cluster, rcut::Float64) = classifyAtoms(cluster.positions, rcut)

function classifyCluster(coordinates::Matrix{Float64}, rcut::Float64)

	tag = classifyAtomsSchebarchov(coordinates, rcut)

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
	lines = readlines(open(filename, "r"))
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

function read_xyz(lines::Vector{SubString{String}})
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

"""
	write_xyz(filename::String, positions::Matrix{Float64})

Takes a path (should end in ".xyz") to a file to write and saves a Cluster type
to that path as a .xyz file.
"""
function write_xyz(filename::String, positions::Matrix{Float64}, formula::Dict{String, Int64}, cell::Float64, tags::Vector{Int64})
	
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
	atoms = read_xyz(filename)
	centreCluster!(atoms)
	write_xyz(filename, atoms, tags)
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

function getStructure(name::String)
	if name == "LJ98_GM"
		xyzCatenation = "98
Lattice=\"15.17457465 0.0 0.0 0.0 15.17457465 0.0 0.0 0.0 15.17457465\" Properties=species:S:1:pos:R:3 pbc=\"F F F\"
Ne       7.54119271       7.47150856      10.14420067
Ne       6.73483947       8.19552962       9.94261875
Ne       8.30690278       8.24649626       9.87930329
Ne       7.52388033       6.48146140       9.64017141
Ne       7.49686125       8.97386463       9.67682288
Ne       6.50428972       6.10999374       9.40141863
Ne       6.72661506       7.20413962       9.42840027
Ne       9.25838047       7.78869936       9.50210427
Ne       8.30414956       7.26050272       9.38476840
Ne       5.91948408       7.92207870       9.23905561
Ne       5.69793995       6.83400652       9.19977837
Ne       8.51889129       9.16269869       9.29468997
Ne       7.49552717       7.98983943       9.17505992
Ne       6.68177915       8.71730543       8.97918611
Ne       7.09100811       9.73487389       8.96038793
Ne       9.26115283       6.80237618       8.99056393
Ne       8.30768452       6.28239486       8.91089130
Ne       9.46940609       8.70074045       8.92014671
Ne       8.50995401       8.17726522       8.79575992
Ne       7.30314322       5.94614354       8.66770715
Ne       7.50064600       7.02031713       8.68489092
Ne       7.70593946       8.89922963       8.59479317
Ne       6.28544305       5.59231478       8.43122280
Ne       6.49525020       6.67644393       8.45848652
Ne       6.69458058       7.74411288       8.48337377
Ne       8.10837857       9.92285577       8.58003653
Ne       5.86834409       8.47280682       8.30104755
Ne       5.68075499       7.40295709       8.26204673
Ne       5.47540311       6.31966132       8.22865421
Ne       9.45974858       7.70908960       8.40644726
Ne       6.27982640       9.47946067       8.24531514
Ne       8.49931524       7.20483843       8.31149443
Ne       8.21498967       5.42877253       8.24119181
Ne       8.72873951       9.09555539       8.21385257
Ne       7.70394215       7.92484083       8.10789104
Ne       7.18740133       5.06947369       8.01857402
Ne       9.94340768       6.74536896       8.11542221
Ne       8.96545453       6.25328660       8.04612994
Ne       6.90172786       8.63939212       7.91211102
Ne       7.89794637       6.41522149       7.88068969
Ne       7.29351692       9.65425961       7.86501899
Ne       7.10001020       7.13182983       7.67522215
Ne       6.89886504       6.06179682       7.64158286
Ne       9.68104193       8.61711212       7.83261254
Ne       8.71027264       8.11658375       7.72972546
Ne       5.02013714       7.01559229       7.47662483
Ne       5.23371437       8.10583704       7.49581054
Ne       6.30038091       7.84976119       7.48127674
Ne       6.09487333       6.78375871       7.44054957
Ne       5.88209317       5.70042291       7.40563361
Ne       7.91874709       8.82732582       7.53186933
Ne       8.30410005       9.85352727       7.48851269
Ne      10.15443751       7.65742679       7.53345469
Ne       9.17077792       7.17082157       7.45809165
Ne       6.02555678       8.89317161       7.31116526
Ne       8.90934240       5.38764874       7.37472410
Ne       8.10295179       7.31716412       7.30023007
Ne       7.81391850       5.54811200       7.22745956
Ne       6.41129041       9.91706922       7.23252597
Ne       9.85428939       7.12881201       6.59136290
Ne       7.30804723       8.03095361       7.10150650
Ne       9.64228660       6.21255068       7.17601266
Ne       6.78587700       5.17709542       6.99768256
Ne       8.94609288       9.04155736       7.15031715
Ne       8.55988516       6.36874587       7.02049141
Ne       7.05022240       9.07500242       6.92800642
Ne       7.49649465       6.52282952       6.86015972
Ne       7.42866705      10.10510096       6.85220722
Ne       5.42496693       6.39918618       6.65737685
Ne       5.64781205       7.49320623       6.68585640
Ne       9.39005832       8.08839725       6.87517812
Ne       6.70501910       7.23353846       6.66226306
Ne       6.49520813       6.16266727       6.61363548
Ne       8.31605483       8.22226050       6.72764922
Ne       5.37731074       8.55927744       6.49162491
Ne       6.17170237       9.32902803       6.30842649
Ne       5.78213786       7.94285765       5.67238617
Ne       7.70364018       7.41809961       6.28890333
Ne       8.99524761       8.80099817       6.08943301
Ne       7.41599404       5.64293412       6.20993017
Ne       8.15745514       6.46986826       5.99211551
Ne       7.19374262       9.51793303       5.92635158
Ne       8.77030786       7.27817901       6.44019811
Ne       6.03910229       6.87931752       5.86565993
Ne       8.50778391       5.49530848       6.35386333
Ne       7.45384012       8.46036567       6.11109680
Ne       9.23882585       6.32072161       6.15049361
Ne       8.07255347       9.27148936       6.54586058
Ne       6.43942118       8.27287429       6.49036972
Ne       6.57838816       8.70976854       5.48542023
Ne       9.47657586       7.83919416       5.81043625
Ne       7.77635439       7.20038848       5.21614562
Ne       8.11180992       9.06470976       5.46928264
Ne       8.38084322       7.99832431       5.66402165
Ne       7.09283752       6.60797932       5.82251794
Ne       8.86397882       7.03475967       5.37156654
Ne       6.84001718       7.65291061       5.66273636
Ne       7.49925062       8.26025112       5.03037398"
	return read_xyz(split(xyzCatenation, "\n"))
	end
end
