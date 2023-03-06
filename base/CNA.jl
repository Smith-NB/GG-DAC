include("Atoms.jl")
"""
	getNeighbourList(atoms::Cluster, rcut::Float64)

Takes a Cluster type and a rcut value (defining bond distance) and returns a table
of all bonds in a cluster.
Return table starts with the number of bonds and then lists the indices of atoms the given
atom is bonded too. After nbond+1 entries, row is padded with 0's.
"""
function getNeighbourList(atoms::Cluster, rcut::Float64)
	r = getDistances(atoms)
	natoms = getNAtoms(atoms)

	#create matrix of bonds between atoms	
	graphbonds = falses(natoms, natoms)
	n1 = 0 
	for i in 1:natoms
		for j in i+1:natoms
			if r[i, j] < rcut
				graphbonds[i, j] = true
				graphbonds[j, i] = true
				n1 += 1
			end
		end
	end

	#create list of bonds
	bondlist = Array{Tuple{Int, Int}}(undef, n1)
	n2 = 0
	for i in 1:natoms
		for j in i+1:natoms
			if graphbonds[i, j]
				n2 += 1
				bondlist[n2] = (i, j)
			end
		end
	end
	return bondlist, graphbonds
end

function DFS(G::BitMatrix, v::Int64, size::Int64, caller::Int64, visited::BitMatrix, commonNeighbours::Array{Int64}, N::Int64)
	if visited[v, 1]
		return size
	end
	visited[v, 1] = true
	for w in 1:N
		if G[v, commonNeighbours[w]] && !visited[commonNeighbours[w], 2] && commonNeighbours[w] != caller
			size = DFS(G, commonNeighbours[w], size+1, v, visited, commonNeighbours, N)
		end
	end
	visited[v, 2] = true
	return size
end



function findLongestChain(graphbonds::BitMatrix, commonNeighbours::Array{Int64}, N::Int64, visited::BitMatrix)
	nl = 0
	#set all common neighbours as unvisited
	for i in 1:N
		visited[commonNeighbours[i], 1] = false
		visited[commonNeighbours[i], 2] = false
	end

	#iterate over all common neighbours and complete DFS if a CN is unvisited. (DFS will be run once per graph component).
	for i in 1:N
		if !visited[commonNeighbours[i], 1]
			size = DFS(graphbonds, commonNeighbours[i], 0, -1, visited, commonNeighbours, N)
			nl = size > nl ? size : nl #set nl index to the component of the largest size.
		end
	end
	return nl
end

function getCNAProfile(atoms::Cluster, rcut::Float64)

	natoms = getNAtoms(atoms)
	bondlist, graphbonds = getNeighbourList(atoms, rcut)
	nbonds = length(bondlist)
	profile = Dict{Tuple{Int64, Int64, Int64}, Int64}()
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)

	#for each bonding pair in cluster
	for i in 1:nbonds
		ncn = 0
		nb = 0
		
		#for each atom in cluster
		for j in 1:natoms
			#check if atom is bonded to both atoms of current pair
			if graphbonds[j, bondlist[i][1]] == true && graphbonds[j, bondlist[i][2]] == true
				ncn += 1
				commonNeighbours[ncn] = j
				#check if bonds exist between current and previosuly discovered common neighbours
				for k in 1:ncn-1
					if graphbonds[commonNeighbours[k], j] == true
						nb += 1
					end
				end
			end
		end
		#find the longest chain of bonds between common neighbours
		nl = findLongestChain(graphbonds, commonNeighbours, ncn, visited)

		#add signature to profile
		sig = (ncn, nb, nl)
		if haskey(profile, sig)
			profile[sig] += 1
		else
			profile[sig] = 1
		end

		#reset oommonNeighbours
		for j in 1:ncn
			commonNeighbours[j] = 0
		end
	end
	return profile
end

function getCNASimilarity(x::Dict{Tuple{Int64, Int64, Int64}, Int64}, y::Dict{Tuple{Int64, Int64, Int64}, Int64})
	intersection = 0
	union = 0
	X = keys(x)	
	for key in X
		a = x[key]
		b = haskey(y, key) ? y[key] : 0
		intersection += a < b ? a : b
		union += a > b ? a : b
	end
	for key in keys(y)
		if key in X
			continue
		end
		a = haskey(x, key) ? x[key] : 0
		b = y[key]
		intersection += a < b ? a : b
		union += a > b ? a : b
	end
	return intersection/union
end
