"""
	isClusterCoherent(clusterCoords::Matrix{Float64}, maxDistance::Number)

Returns true if a cluster is coherent, else false. Cohenerncy is that all atoms form one cluster - if one or more atoms
	cannot be reached by any other atom in a graph via bonds of size maxDistance, the cluster is incoherent.
"""
function isClusterCoherent(clusterCoords::Matrix{Float64}, maxDistance::Number)

	N =	getNAtoms(clusterCoords)
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
function generateRandomSeed(formula::Dict{String, Int64}, boxLength::Number, vacuumAdd::Number, returnCoordsOnly::Bool=false, coherencyDistance::Float64=1.5)
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
		if isClusterCoherent(clusterCoords, coherencyDistance)
			break
		else
			clusterCoords = zeros(Float64, N, 3)
			println("Cluster non-coherent. Trying again. $(Threads.threadid())")
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
	
	if coords[75, 1] > 100 || coords[75, 1] < -100 || coords[75, 2] > 100 || coords[75, 2] < -100 || coords[75, 3] > 100 || coords[75, 3] < -100
		println("COORDALERT1")
	end

	# get radius to place cluster on to.
	radius = getInclusionRadiusOfCluster(coords)
	if radius == Inf || radius == -Inf
		println("INFALERT_RADIUS")
	end
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
		if Inf in newCoords[atomsToMove[i], :] || -Inf in newCoords[atomsToMove[i], :]
			println("INFALERT3")
			println(centreOfMass)
			println("$radius $sphericalCoords")
			println(sphericalToCartesian(radius, sphericalCoords[i, 1], sphericalCoords[i, 2]))
		end
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
	perturbClusterSurface!(atoms::Cluster, nAtomsToMove::Number, rCut::Float64)

Sends atom positions to perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)
"""
function perturbClusterSurface!(atoms::Cluster, nAtomsToMove::Number, rCut::Float64)
	# get radius to place cluster on to.
	radius = getInclusionRadiusOfCluster(atoms.positions)
	centreOfMass = getCentreOfCluster(atoms.positions)
	N = getNAtoms(atoms.positions)
	if nAtomsToMove < 1
		nAtomsToMove = trunc(Int64, round(nAtomsToMove*N, digits=2))
	end

	# determine lowest coordination number atoms
	r = atoms.distances
	
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
			index = minCoordAtoms[rand(0x0d:length(minCoordAtoms))]
			while index in atomsToMove
				index = minCoordAtoms[rand(0x0d:length(minCoordAtoms))]
			end
			atomsToMove[i] = index
		end
	end

	sphericalCoords = rand(nAtomsToMove, 2)*360
	
	for i in 1:nAtomsToMove
		atoms.positions[atomsToMove[i], :] = sphericalToCartesian(radius, sphericalCoords[i, 1], sphericalCoords[i, 2])
		atoms.positions[atomsToMove[i], 1] += centreOfMass[1]
		atoms.positions[atomsToMove[i], 2] += centreOfMass[2]
		atoms.positions[atomsToMove[i], 3] += centreOfMass[3]
	end

	return nothing
end

"""
	perturbCluster(coords::Matrix{Float64}, dr::Float64)

Returns positions of atoms after moving each atom in each coordinate direction by ±`dr`
"""
function perturbCluster(coords::Matrix{Float64}, dr::Float64)
	n = getNAtoms(coords)
	#return coords + rand(-1.0*dr:0.00001*dr:1.0*dr, n, 3)
	r = rand(n, 3)
	rCopy = copy(r)
	for i in 1:n
		r[i, 1] = (2*r[i, 1] - 1)*dr + coords[i, 1]
		r[i, 2] = (2*r[i, 2] - 1)*dr + coords[i, 2]
		r[i, 3] = (2*r[i, 3] - 1)*dr + coords[i, 3]
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
	perturbCluster!(coords::Matrix{Float64}, dr::Float64)

Sends atom positions to perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Number, rCut::Float64)
"""
function perturbCluster!(atoms::Cluster, dr::Float64)
	n = getNAtoms(atoms)
	#return coords + rand(-1.0*dr:0.00001*dr:1.0*dr, n, 3)
	for i in 1:n
		atoms.positions[n, 1] += rand()*dr*2-1
		atoms.positions[n, 2] += rand()*dr*2-1
		atoms.positions[n, 3] += rand()*dr*2-1
	end
	return nothing
end

function geometricCentreDisplacement!(atoms::Cluster, alphaMin::Float64, alphaMax::Float64, w::Float64)
	CoM = getCentreOfCluster(atoms)
	N = getNAtoms(atoms)
	R = Vector{Float64}(undef, N)
	rho = Vector{Float64}(undef, N)
	alphaDiff = alphaMax - alphaMin

	for i in 1:N
		d = atoms.positions[i, 1] - CoM[1]
		d += atoms.positions[i, 2] - CoM[2]
		d += atoms.positions[i, 3] - CoM[3]
		R[i] = d
	end
	
	Rmax = maximum(R)
	for i in 1:N
		rho[i] = alphaDiff * (R[i]/Rmax)^w + alphaMin
	end

	polar = rand(N) * 2*π		# θ
	azimuth = rand(N) * 2*π 	# φ
	for i in 1:N
		atoms.positions[i, 1] += rho[i] * sin(polar[i]) * cos(azimuth[i])
		atoms.positions[i, 2] += rho[i] * sin(polar[i]) * sin(azimuth[i])
		atoms.positions[i, 3] += rho[i] * cos(polar[i])
	end

end

function geometricCentreDisplacement(coords::Matrix{Float64}, alphaMin::Float64, alphaMax::Float64, w::Float64)
	CoM = getCentreOfCluster(coords)
	N = getNAtoms(coords)
	R = Vector{Float64}(undef, N)
	rho = Vector{Float64}(undef, N)
	alphaDiff = alphaMax - alphaMin

	for i in 1:N
		d = coords[i, 1] - CoM[1]
		d += coords[i, 2] - CoM[2]
		d += coords[i, 3] - CoM[3]
		R[i] = d
	end
	
	Rmax = maximum(R)
	for i in 1:N
		rho[i] = alphaDiff * (R[i]/Rmax)^w + alphaMin
	end

	polar = rand(N) * 2*π		# θ
	azimuth = rand(N) * 2*π 	# φ
	newPositions = Matrix{Float64}(undef, N, 3)
	for i in 1:N
		newPositions[i, 1] = coords[i, 1] + rho[i] * sin(polar[i]) * cos(azimuth[i])
		newPositions[i, 2] = coords[i, 2] + rho[i] * sin(polar[i]) * sin(azimuth[i])
		newPositions[i, 3] = coords[i, 3] + rho[i] * cos(polar[i])
	end

	return newPositions
end

geometricCentreDisplacement(atoms::Cluster, alphaMin::Float64, alphaMax::Float64, w::Float64) = geometricCentreDisplacement(atoms.positions, alphaMin, alphaMax, w)


"""	
	getSeedFromPool(pool::Vector{Matrix{Float64}}, n::Int64) 
	
Given a `pool` of cluster coordinates of size `n`, randomly return one the coordinates of one
	cluster in the pool.

"""
function getSeedFromPool(pool::Vector{Matrix{Float64}}, n::Int64)
	r = rand(0x0d:n)
	return pool[r]
end