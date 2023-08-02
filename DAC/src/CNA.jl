
const CNAProfile = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}

#struct CNAProfile <: Dict{Tuple{UInt8, UInt8, UInt8}, UInt16} end

"""
	getNeighbourList(atoms::Cluster, rcut::Float64)

Takes a Cluster type and a rcut value (defining bond distance) and returns a table
of all bonds in a cluster.
Return table starts with the number of bonds and then lists the indices of atoms the given
atom is bonded too. After nbond+1 entries, row is padded with 0's.
"""
function getNeighbourList(coordinates::Matrix{Float64}, rcut::Float64)
	r = getDistances(coordinates)
	natoms = getNAtoms(coordinates)

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

getNeighbourList(atoms::Cluster, rcut::Float64) = getNeighbourList(atoms.positions, rcut)


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

function getNormalCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
	natoms = getNAtoms(coordinates)
	bondlist, graphbonds = getNeighbourList(coordinates, rcut)
	nbonds = length(bondlist)
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)

	normalCNA = Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, natoms)
	for i in 1:natoms
		normalCNA[i] = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}() 
	end

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
		if haskey(normalCNA[bondlist[i][1]], sig)
			normalCNA[bondlist[i][1]][sig] += 1
		else
			normalCNA[bondlist[i][1]][sig] = 1
		end

		if haskey(normalCNA[bondlist[i][2]], sig)
			normalCNA[bondlist[i][2]][sig] += 1
		else
			normalCNA[bondlist[i][2]][sig] = 1
		end

		#reset oommonNeighbours
		for j in 1:ncn
			commonNeighbours[j] = 0
		end
	end

	return normalCNA

end

getNormalCNAProfile(atoms::Cluster, rcut::Float64) = getNormalCNAProfile(atoms.positions, rcut)


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


"""
	getCNASimilarity(x::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}, y::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})

Takes two cluster CNA profiles (formatted as dictionaries) and calculates
the similarity between them.
"""
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

function CNAToString(cna::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})
	s = ""
	for sig in keys(cna)
		s *= "($(sig[1]), $(sig[2]), $(sig[3])): $(cna[sig]) "
	end
	return s
end

function CNAToString(cna::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}})
	s = ""
	for i in 1:length(cna)
		s *= "($(cna[i].first[1]), $(cna[i].first[2]), $(cna[i].first[3])): $(cna[i].second) "
	end
	return s
end

function CNAToLogString(cna::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}})
	# create the string to log from the given CNA and cluster ID
	s = ""
	for pair in cna
		s *= string(pair.first[1]) * "," * string(pair.first[2]) * "," * string(pair.first[3]) * ":" * string(pair.second) * ";"
	end
	return s
end
"""
	stringToCNA(s::String)

Takes a string formatted as: "ncn,nb,nl:freq;ncn,nb,nl:freq;...ncn,nb,nl:freq;"
and converts it to a CNA profile of type 
Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}} 
where:
	ncn,nb,nl is the CNA signagure, e.g. (5, 5, 5) = "5,5,5"
	freq is the frequency 
"""
function stringToCNA(s::String)
	
	freqCNAPair = [split(x, ':') for x in split(s, ';')]
	n = length(freqCNAPair[length(freqCNAPair)]) == 1 ? length(freqCNAPair) -1 : length(freqCNAPair)
	CNA = Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, n)
	
	# For all CNA signature & frequency pairs
	for i in 1:n
		cna::Tuple{UInt8, UInt8, UInt8} = Tuple(parse(UInt8, String(x)) for x in split(freqCNAPair[i][1], ','))
		freq::UInt16 = parse(UInt16, String(freqCNAPair[i][2]))
		CNA[i] = Pair(cna, freq)
	end

	return CNA
end

function getCNA(name::String)
	if name == "LJ38_GM"
		return "4,2,1:60;3,1,1:48;2,1,1:24;2,0,0:12;"
	elseif name == "LJ55_GM"
		return "5,5,5:24;4,2,2:90;3,2,2:60;3,1,1:60;"
	elseif  name == "LJ75_GM"
		return "5,5,5:4;4,2,2:65;4,2,1:90;3,2,2:20;3,1,1:90;3,0,0:10;2,1,1:20;2,0,0:20;"
	elseif  name == "LJ75_Ih"
		return "5,5,5:42;4,3,3:14;4,2,2:133;3,2,2:55;3,1,1:55;3,0,0:2;2,1,1:7;2,0,0:20;"
	elseif  name == "LJ98_GM"
		return "5,5,5:18;4,3,3:24;4,2,2:156;4,2,1:66;3,2,2:12;3,1,1:72;3,0,0:12;2,1,1:36;2,0,0:36;"
	elseif  name == "LJ104_GM"
		return "5,5,5:5;4,2,2:97;4,2,1:158;3,2,2:22;3,1,1:106;3,0,0:12;2,1,1:41;2,0,0:16;1,0,0:2;"
	elseif  name == "LJ98_FCC1"
		return "4,2,2:42;4,2,1:204;3,2,2:9;3,1,1:91;3,0,0:2;2,1,1:59;2,0,0:15;1,0,0:2;"
	elseif  name == "LJ98_FCC2"
		return "5,5,5:3;4,2,2:70;4,2,1:173;3,2,2:17;3,1,1:87;3,0,0:2;2,1,1:50;2,0,0:23;1,0,0:1;"
	elseif name == "LJ98_Dh"
		return "5,5,5:4;4,2,2:85;4,2,1:143;3,2,2:15;3,1,1:101;3,0,0:18;2,1,1:35;2,0,0:27;"
	end
end