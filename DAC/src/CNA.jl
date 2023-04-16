
#struct CNAProfile <: Dict{Tuple{UInt8, UInt8, UInt8}, UInt16} end

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
	"""
	for i in 1:natoms
		println(r[i, :])
	end
	"""
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


"""
	DFS(G::BitMatrix, v::Int64, size::Int64, caller::Int64, visited::BitMatrix, commonNeighbours::Array{Int64}, N::UInt8)

Performs a depth first search on a graph subset starting from common neighbour v. Marks all atoms encountered as visited.
"""
function DFS(G::BitMatrix, v::Int64, size::Int64, caller::Int64, visited::BitMatrix, commonNeighbours::Array{Int64}, N::UInt8)
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


"""
	findLongestChain(graphbonds::BitMatrix, commonNeighbours::Array{Int64}, N::UInt8, visited::BitMatrix)

Takes a BitMatrix specifiying the bonding network of a cluster, and a set of N common neighbours of an arbitrary bonding pair.
The longest chain of bonds, nl, between these common neighbours will then be returned.

"""
function findLongestChain(graphbonds::BitMatrix, commonNeighbours::Array{Int64}, N::UInt8, visited::BitMatrix)
	nl::UInt8 = 0
	#set all common neighbours as unvisited
	for i in 1:N
		visited[commonNeighbours[i], 1] = false
		visited[commonNeighbours[i], 2] = false
	end

	#iterate over all common neighbours and complete DFS if a common neighbour is unvisited. (DFS will be run once per graph component).
	for i in 1:N
		if !visited[commonNeighbours[i], 1]
			size = DFS(graphbonds, commonNeighbours[i], 0, -1, visited, commonNeighbours, N)
			nl = size > nl ? size : nl #set nl index to the component of the largest size.
		end
	end
	return nl
end

function getCNAProfileAsDict(atoms::Cluster, rcut::Float64)
	returnDict = true

	natoms = getNAtoms(atoms)
	bondlist, graphbonds = getNeighbourList(atoms, rcut)
	#println("bondlist ", bondlist)
	nbonds = length(bondlist)
	profile = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}()
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)
	nsigs = 0

	#for each bonding pair in cluster
	for i in 1:nbonds
		ncn::UInt8 = 0
		nb::UInt8 = 0

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
		nl::UInt8 = findLongestChain(graphbonds, commonNeighbours, ncn, visited)

		#add signature to profile
		sig = (ncn, nb, nl)
		if haskey(profile, sig)
			profile[sig] += 1
		else
			nsigs += 1
			profile[sig] = 1
		end

		#reset oommonNeighbours
		for j in 1:ncn
			commonNeighbours[j] = 0
		end
	end
	
	if returnDict
		return profile
	end

	storedCNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, length(profile))
	CNAkeys = sort(collect(keys(profile)), rev=true)
	for i in 1:length(CNAkeys)
		storedCNA[i] = Pair(CNAkeys[i], profile[CNAkeys[i]])
		#println(storedCNA[i])
	end

	return storedCNA
end

function getCNAProfile(atoms::Cluster, rcut::Float64)
	natoms = getNAtoms(atoms)
	bondlist, graphbonds = getNeighbourList(atoms, rcut)
	#println("bondlist ", bondlist)
	nbonds = length(bondlist)
	profile = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}()
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)
	nsigs = 0

	#for each bonding pair in cluster
	for i in 1:nbonds
		ncn::UInt8 = 0
		nb::UInt8 = 0

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
		nl::UInt8 = findLongestChain(graphbonds, commonNeighbours, ncn, visited)

		#add signature to profile
		sig = (ncn, nb, nl)
		if haskey(profile, sig)
			profile[sig] += 1
		else
			nsigs += 1
			profile[sig] = 1
		end

		#reset oommonNeighbours
		for j in 1:ncn
			commonNeighbours[j] = 0
		end
	end

	storedCNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, length(profile))
	CNAkeys = sort(collect(keys(profile)), rev=true)
	for i in 1:length(CNAkeys)
		storedCNA[i] = Pair(CNAkeys[i], profile[CNAkeys[i]])
		#println(storedCNA[i])
	end

	return storedCNA
end

function getCNASimilarity(x::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}, y::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}})
	intersection = 0
	union = 0
	Nx = length(x)
	Ny = length(y)


	for i in 1:Nx
		sig = x[i].first
		y_index = binarySearch(y, Ny, sig)

		a = x[i].second								# get frequency of sig in cluster x
		b = y_index < 0 ? 0 : y[y_index].second		# get frequency of sig in cluster y, or 0 if sig is not present.
		intersection += a < b ? a : b 				# add the smaller of a or b to the intersection.
		union += a > b ? a : b 						# add the larger of a or b to the union.

	end

	for i in 1:Ny
		sig = y[i].first
		x_index = binarySearch(x, Nx, sig)
		union += x_index == -1 ? y[i].second : 0
	end

	return intersection / union



	# The following is left for inspiration, and involved a queue implementation of the above.
	# This way, only unchecked sigs in y needed to be iterated through, those sigs being the
	# ones remaining in the queue. However this requried memory allocations, which were not cost effective.
	#=
	
	sigIndexChecked_x = collect(1:Nx)
	sigIndexChecked_y = collect(1:Ny)

	# parse through all signatures of cluster x
	while length(sigIndexChecked_x) != 0
		i = pop!(sigIndexChecked_x) 	# get index of signature in x
		sig = x[i].first				# get signature

		y_index = binarySearch(y, Ny, sig) 		# look for selected signature in other cluster.
		removeall!(sigIndexChecked_y, y_index)	# mark signature as checked.

		a = x[i].second								# get frequency of sig in cluster x
		b = y_index < 0 ? 0 : y[y_index].second		# get frequency of sig in cluster y, or 0 if sig is not present.
		intersection += a < b ? a : b 				# add the smaller of a or b to the intersection.
		union += a > b ? a : b 						# add the larger of a or b to the union.
	end

	# parse through remaining, unchecked signatues of cluster y.
	# only need to add sig frequencies to union, as these signatures can't be present in cluster x,
	# i.e. intersection of all remaining signatures in cluster y is 0.
	while length(sigIndexChecked_y) != 0
		i = pop!(sigIndexChecked_y)
		union += y[i].second
	end

	return intersection / union
	=#
end

function getCNASimilarity(x::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}, y::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})
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

function stringToCNA(s::String)
	#5,5,5:16;5,4,4:8;4,3,3:26;4,2,2:15;4,2,1:7;3,2,2:21;3,1,1:16;3,0,0:1;2,1,1:20;2,0,0:14;1,0,0:1;
	#Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}
	freqCNAPair = [split(x, ':') for x in split(s, ';')]
	CNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, length(freqCNAPair))
	for i in 1:length(freqCNAPair)
		cna::Tuple{UInt8, UInt8, UInt8} = (parse(UInt8, x) for x in freqCNAPair[i][1])
		freq::UInt16 = freqCNAPair[i][2]
		CNA[i] = Pair(cna, freq)
	end
end