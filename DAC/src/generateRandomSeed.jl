"""
	isClusterCoherent(clusterCoords::Matrix{Float64}, maxDistance::Number)

Returns true if a cluster is coherent, else false. Cohenerncy is that all atoms form one cluster - if one or more atoms
	cannot be reached by any other atom in a graph via bonds of size maxDistance, the cluster is incoherent.
"""
function isClusterCoherent(clusterCoords::Matrix{Float64}, maxDistance::Number)

	N = trunc(Int, length(clusterCoords)/3)
	neighboursList = getNeighboursList(clusterCoords, maxDistance)
	atomsNotReached = [i for i in 1:N]
	clusterPaths = [1]
	while length(clusterPaths) != 0
		atomToExplore = pop!(clusterPaths)
		index = findfirst([i==atomToExplore for i in atomsNotReached])
		if index != nothing
			deleteat!(atomsNotReached, index)
		end
		append!(clusterPaths, neighboursList[atomToExplore])
		while length(neighboursList[atomToExplore]) != 0
			atomBondedTo = neighboursList[atomToExplore][1]
			removefirst!(neighboursList[atomBondedTo], atomToExplore)
			removefirst!(neighboursList[atomToExplore], atomBondedTo)
		end
	end

	return length(atomsNotReached) == 0
end

"""
	generateRandomSeed(formula::Dict{String, Int64}, boxLength::Number, vacuumAdd::Number, returnCoordsOnly::Bool=false)

Returns a DC.Cluster or Matrix{Float64} (`returnCoordsOnly` dependant) type with randomly generated positions of atoms inside
	a box of size `boxLength` and with a cell of size `boxLength` + `vacuumAdd`.
"""
function generateRandomSeed(formula::Dict{String, Int64}, boxLength::Number, vacuumAdd::Number, returnCoordsOnly::Bool=false)
	# get number of atoms
	N = 0
	for key in keys(formula)
		N += formula[key]
	end


	clusterCoords = zeros(Float64, N, 3)

	while true
		# Generate N atomic coordinates.
		for i in 1:N
			# Generate xyz coords for new atom and check if it is too close to other atoms.
			atomCoords = rand(Float64, 3)
			while true
				atomClash = false
				clusterCoords[i, 1:3] = rand(Float64, 3) * boxLength

				for j in 1:i-1
					if abs(clusterCoords[i, 1] - clusterCoords[j, 1]) < 0.0001 && abs(clusterCoords[i, 2] - clusterCoords[j, 2]) < 0.0001 && abs(clusterCoords[i, 3] - clusterCoords[j, 3]) < 0.0001
						atomClash = true
					end
				end

				if !atomClash
					break
				end
			end
		end

		# Check that the cluster is coherent. Restart if it isn't. 
		if isClusterCoherent(clusterCoords, 1.5)
			break
		else
			clusterCoords = zeros(Float64, N, 3)
			println("Cluster non-coherent. Trying again.")
		end
	end
	
	if returnCoordsOnly
		return clusterCoords
	else	
		cellLength = boxLength + vacuumAdd
		cell = zeros(Float64, 3, 3)
		cell[1] = cellLength; cell[5] = cellLength; cell[9] = cellLength
		return Cluster(formula, clusterCoords, cell)
	end

end

"""
	perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)

Returns new postions of a cluster after moving `nAtomsToMove` atoms on the surface.
	`nAtomsToMove` can be a number (number of atoms) or float < 1 (% of total atoms).
	Only the least coordinated atoms are selected. If the number of lowest coordination atoms 
	is less than `nAtomsToMove` then only this lower number of atoms are displaced.
	Atoms are moved to a random location on the cluster surface.
"""
function perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)
	
	# get radius to place cluster on to.
	radius = getInclusionRadiusOfCluster(coords)
	centreOfMass = getCentreOfCluster(coords)
	N = getNAtoms(coords)
	if nAtomsToMove < 1
		nAtomsToMove = trunc(Int64, round(nAtomsToMove*N, digits=2))
	end

	# determine lowest coordination number atoms
	r = getDistances(coords)
	
	coordinationNumbers  = zeros(Int64, N)
	for i in 1:N
		for j in i+1:N
			if r[j, i] > 0 && r[j, i] <= rCut
				coordinationNumbers[i] += 1
				coordinationNumbers[j] += 1
			end
		end
	end
	minCoord::Int64 = minimum(coordinationNumbers)
	minCoordAtoms = findall(i->i==minCoord, coordinationNumbers)

	atomsToMove = nothing
	# if the number of low coord atoms is less than nAtomsToMove, only move the fewer amount
	if length(minCoordAtoms) <= nAtomsToMove 
		nAtomsToMove = length(minCoordAtoms)
		atomsToMove = minCoordAtoms
	else # otherwise randomly select which atoms to move
		atomsToMove = Vector{Int64}(undef, nAtomsToMove)
		for i in 1:nAtomsToMove
			index = minCoordAtoms[rand(1:length(minCoordAtoms))]
			while index in atomsToMove
				index = minCoordAtoms[rand(1:length(minCoordAtoms))]
			end
			atomsToMove[i] = index
		end
	end

	sphericalCoords = rand(nAtomsToMove, 2)*360
	newCoords = copy(coords)
	for i in 1:nAtomsToMove
		newCoords[atomsToMove[i], :] = sphericalToCartesian(radius, sphericalCoords[i, 1], sphericalCoords[i, 2])
		newCoords[atomsToMove[i], 1] += centreOfMass[1]
		newCoords[atomsToMove[i], 2] += centreOfMass[2]
		newCoords[atomsToMove[i], 3] += centreOfMass[3]
	end

	return newCoords

end

"""
	perturbClusterSurface(atoms::Cluster, nAtomsToMove::Number, rCut::Float64)

Sends atom positions to perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)
"""
function perturbClusterSurface(atoms::Cluster, nAtomsToMove::Number, rCut::Float64)
	return perturbClusterSurface(atoms.positions, nAtomsToMove, rCut)
end

"""
	perturbCluster(coords::Matrix{Float64}, dr::Float64)

Returns positions of atoms after moving each atom in each coordinate direction by ±`dr`
"""
function perturbCluster(coords::Matrix{Float64}, dr::Float64)
	n = getNAtoms(coords)
	#return coords + rand(-1.0*dr:0.00001*dr:1.0*dr, n, 3)
	r = rand(n, 3)
	for i in 1:n
		r[n, 1] = 2*r[n, 1] + coords[n, 1] - 1
		r[n, 2] = 2*r[n, 2] + coords[n, 2] - 1
		r[n, 3] = 2*r[n, 3] + coords[n, 3] - 1
	end
	return r
end

"""
	perturbCluster(coords::Matrix{Float64}, dr::Float64)

Sends atom positions to perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)
"""
function perturbCluster(atoms::Cluster, dr::Float64)
	return perturbCluster(atoms.positions, dr)
end

"""	
	getSeedFromPool(pool::Vector{Matrix{Float64}}, n::Int64) 
	
Given a `pool` of cluster coordinates of size `n`, randomly return one the coordinates of one
	cluster in the pool.

"""
function getSeedFromPool(pool::Vector{Matrix{Float64}}, n::Int64)
	r = rand(1:n)
	return pool[r]
end