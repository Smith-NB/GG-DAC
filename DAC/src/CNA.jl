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

#getNeighbourList(atoms::Cluster, rcut::Float64) = getNeighbourList(atoms.positions, rcut)

"""
	getNeighbourList(atoms::Cluster, rcut::Float64)

Takes a Cluster type and a rcut value (defining bond distance) and returns a table
of all bonds in a cluster.
Return table starts with the number of bonds and then lists the indices of atoms the given
atom is bonded too. After nbond+1 entries, row is padded with 0's.
"""
function getNeighbourList(atoms::Cluster, rcut::Float64)
	r = atoms.distances
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

"""
	getCNAProfileAsDict(atoms::Cluster, rcut::Float64)

Calculates the total CNA profile of atoms and returns it as a Dict type (hash table).
"""
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

	storedCNA = CNAProfile(undef, length(profile))
	CNAkeys = sort(collect(keys(profile)), rev=true)
	for i in 1:length(CNAkeys)
		storedCNA[i] = Pair(CNAkeys[i], profile[CNAkeys[i]])
		#println(storedCNA[i])
	end

	return storedCNA
end


"""
	getCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)

Calculates the total CNA profile for the cluster described by coordinates
and returns it as a sorted Vector.
"""
function getCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
	natoms = getNAtoms(coordinates)
	bondlist, graphbonds = getNeighbourList(coordinates, rcut)
	nbonds = length(bondlist)
	#profile = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}()
	storedCNA = CNAProfile(undef, 0)
	storedCNA_freq = Vector{UInt16}(undef, 0)
	n_storedCNA = 0
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
		sig::Tuple{UInt8, UInt8, UInt8} = (ncn, nb, nl)
		index::Int64 = binarySearch(storedCNA, n_storedCNA, sig)
		if index < 0
			insert!(storedCNA, -index, Pair(sig, 1))
			insert!(storedCNA_freq, -index, 1)
			n_storedCNA += 1
		else
			storedCNA_freq[index] += 1
		end
		
	end

	
	for i in 1:length(storedCNA)
		storedCNA[i] = Pair(storedCNA[i].first, storedCNA_freq[i])
	end
	
	return storedCNA
end


"""
	getCNAProfile(atoms::Cluster, rcut::Float64)

Calculates the total CNA profile for the cluster described by coordinates
and returns it as a sorted Vector.
"""
getCNAProfile(atoms::Cluster, rcut::Float64) = getCNAProfile(atoms.positions, rcut)

"""
	getCNAProfile(atoms::Cluster, rcut::Float64)

Calculates the total and normal CNA profiles for the cluster described by `coordinates`
and returns it as a sorted `Vector` and a `Vector` of `Dict`'s, respectively
"""
function getTotalAndNormalCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
	natoms = getNAtoms(coordinates)
	bondlist, graphbonds = getNeighbourList(coordinates, rcut)
	nbonds = length(bondlist)
	#profile = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}()
	totalCNA = CNAProfile(undef, 0)
	totalCNA_freq = Vector{UInt16}(undef, 0)
	n_totalCNA = 0
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)
	nsigs = 0

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
		sig::Tuple{UInt8, UInt8, UInt8} = (ncn, nb, nl)
		index::Int64 = binarySearch(totalCNA, n_totalCNA, sig)
		if index < 0
			insert!(totalCNA, -index, Pair(sig, 1))
			insert!(totalCNA_freq, -index, 1)
			n_totalCNA += 1
		else
			totalCNA_freq[index] += 1
		end

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
		
	end

	
	for i in 1:length(totalCNA)
		totalCNA[i] = Pair(totalCNA[i].first, totalCNA_freq[i])
	end
	
	return totalCNA, normalCNA
end


"""
	getTotalAndNormalCNAProfile(atoms::Cluster, rcut::Float64)

Calculates the total and normal CNA profiles for the cluster `atoms`
and returns it as a sorted `Vector` and a `Vector` of `Dict`'s, respectively
"""
getTotalAndNormalCNAProfile(atoms::Cluster, rcut::Float64) = getTotalAndNormalCNAProfile(atoms.positions, rcut)

"""
	getNormalCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)

Calculates the normal CNA profiles for the cluster described by `coordinates`
and returns it as a `Vector` of `Dict`'s, respectively.
`normalCNA[1]` is the atom level CNA profile of the atom at `coordinates[1, :]`
"""
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

"""
	getNormalCNAProfile(atoms::Cluster, rcut::Float64)

Calculates the normal CNA profiles for the cluster described by `coordinates`
and returns it as a `Vector` of `Dict`'s, respectively.
`normalCNA[1]` is the atom level CNA profile of the atom at 
`atom.positions[1, :]`
"""
getNormalCNAProfile(atoms::Cluster, rcut::Float64) = getNormalCNAProfile(atoms.positions, rcut)


function getNormalCNAProfileAsVector(coordinates::Matrix{Float64}, rcut::Float64)
	natoms = getNAtoms(coordinates)
	bondlist, graphbonds = getNeighbourList(coordinates, rcut)
	nbonds = length(bondlist)
	commonNeighbours = zeros(Int64, natoms)
	visited = trues(natoms, 2)

	normalCNA = Vector{CNAProfile}(undef, natoms)
	normalCNA_freq = Vector{Vector{UInt16}}(undef, natoms)
	n_normalCNA = zeros(Int64, natoms)
	
	for i in 1:natoms
		normalCNA[i] = CNAProfile(undef, 0) 
		normalCNA_freq[i] = Vector{UInt16}(undef, 0) 
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
		index1::Int64 = binarySearch(normalCNA[bondlist[i][1]], n_normalCNA[bondlist[i][1]], sig)
		index2::Int64 = binarySearch(normalCNA[bondlist[i][2]], n_normalCNA[bondlist[i][2]], sig)
		if index1 < 0
			insert!(normalCNA[bondlist[i][1]], -index1, Pair(sig, 1))
			insert!(normalCNA_freq[bondlist[i][1]], -index1, 1)
			n_normalCNA[bondlist[i][1]] += 1
		else
			normalCNA_freq[bondlist[i][1]][index1] += 1
		end

		if index2 < 0
			insert!(normalCNA[bondlist[i][2]], -index2, Pair(sig, 1))
			insert!(normalCNA_freq[bondlist[i][2]], -index2, 1)
			n_normalCNA[bondlist[i][2]] += 1
		else
			normalCNA_freq[bondlist[i][2]][index2] += 1
		end


		#reset oommonNeighbours
		for j in 1:ncn
			commonNeighbours[j] = 0
		end
	end

	for i in 1:natoms
		for j in 1:length(normalCNA_freq[i])
			normalCNA[i][j] = Pair(normalCNA[i][j].first, normalCNA_freq[i][j]) # line causes no additional allocs.
		end
	end

	return normalCNA
end

getNormalCNAProfile(atoms::Cluster, rcut::Float64) = getNormalCNAProfile(atoms.positions, rcut)

"""
	getCNASimilarity(x::CNAProfile, y::CNAProfile)

Returns as a Float (0-1) the similarity between two total CNA profiles, `x` and `y`
according to the Jaccard similarity index.
"""
function getCNASimilarity(x::CNAProfile, y::CNAProfile)
	intersection = 0
	union = 0
	Nx = length(x)
	Ny = length(y)

	if Nx < Ny
		temp = x
		x = y
		y = temp
		Nx = length(x)
		Ny = length(y)
	end

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


"""
	getClasses(s::String)

Returns a `Vector` of `Dict`'s of atom level CNA profiles to define `n` classes,
where `n` is the number of lines in the input `String` `s`.
Format example: '421 421 421\n555 555 555' is two classes each consisting of three occurences of
two signatures, (4, 2, 1) and (5, 5, 5), respectively.
"""
function getClasses(s::String)
	lines = split(s, '\n')
    # vector to hold class defining sets of CNA signatures
    classes = Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}(undef, length(lines))
    for n in 1:length(lines)
        sigs = split(lines[n])
        uSigs = unique(sigs) # get set of all unique signatures for class
        cSigs = [count(x->x==u, sigs) for u in uSigs] # get corresponding frequency of each signature
        classes[n] = Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}() # instantiate entry in classes
        for i in 1:length(uSigs)
            # set dictionary keys and values for each unique signature (parse converts String/chars to given number datatype).
            classes[n][(parse(UInt8, uSigs[i][1]), parse(UInt8, uSigs[i][2]), parse(UInt8, uSigs[i][3]))] = cSigs[i]
        end
    end

    return classes
end


"""
	getClasses()

Returns a `Vector` of `Dict`'s of atom level CNA profiles to define the 63 classes
used by Roncaglia & Ferrando.
"""
function getClasses()

    # specify classes
    s = """
    444 444 444 444 444 444 666 666 666 666 666 666 666 666
    433 433 433 433 433 433 555 555 555 555 555 555 666 666
    444 444 444 444 444 544 544 544 544 666 666 666 666
    433 433 433 433 433 433 555 555 555 555 555 555 666
    421 421 421 421 421 421 421 421 421 421 421 421
    421 421 421 421 421 421 422 422 422 422 422 422
    555 555 555 555 555 555 555 555 555 555 555 555
    422 422 422 422 422 422 422 422 422 422 555 555
    422 422 422 422 422 422 433 433 544 544 555 555
    421 421 421 421 422 422 422 422 433 433 544 544
    421 421 422 422 422 422 433 433 544 544 555 555
    300 300 311 311 422 422 433 433 544 544 555 555
    311 311 311 311 421 421 421 421 421 421 421
    311 311 311 311 421 421 421 422 422 422 422
    300 300 300 300 300 422 422 422 422 422 555
    300 311 311 322 421 421 421 421 422 422 422
    311 311 311 311 421 421 421 421 421 422 422
    211 300 300 311 311 421 421 422 422 422 422
    200 300 300 300 300 322 422 422 422 422 555
    433 433 433 433 433 555 555 555 555 555 555
    211 311 311 311 311 421 421 421 421 421
    211 311 311 311 311 421 421 421 422 422
    322 322 433 433 433 433 444 444 666 666
    300 300 311 311 311 311 421 421 422 422
    311 311 311 311 422 422 422 422 433 433
    311 311 311 311 433 433 433 433 544 544
    200 200 300 300 322 422 422 422 422 555
    200 300 300 311 422 422 422 433 433 555
    200 300 300 300 322 322 422 422 422 555
    311 311 311 311 311 311 421 421 421
    211 311 311 322 322 421 421 422 422
    211 211 211 211 444 544 544 544 544
    311 311 311 311 322 421 433 433 544
    311 311 322 322 433 433 433 544 555
    322 322 322 433 433 433 555 555 555
    200 200 200 200 200 200 211 211 211
    211 211 211 211 421 421 421 421
    200 200 211 311 311 421 422 422
    311 311 311 311 322 322 422 422
    200 200 322 322 322 422 422 555
    311 311 311 311 322 322 422 422
    322 322 322 322 433 433 444 666
    211 211 211 433 433 444 544 544
    322 322 322 322 433 433 555 555
    200 200 311 311 311 311 421
    211 211 211 311 311 421 421
    200 211 211 311 322 421 422
    200 211 211 311 322 421 422
    200 200 300 311 311 322 422
    200 200 211 211 422 422 433
    200 200 300 311 311 322 422
    322 322 322 322 322 322 666
    211 211 322 433 433 444 544
    211 211 311 311 311 311
    100 100 211 211 422 422
    200 211 211 311 311 421
    322 322 322 322 322 555
    322 322 322 322 444
    100 211 211 322 422
    200 200 211 311 311
    200 200 200 200
    211 211 322 322
    211 211 211
    """

    return getClasses(s)
end

"""
	getAtomClasses(nCNA::Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}, classes::Vector{Dict{Tuple{UInt8, UInt8, UInt8}, UInt16}}, N::Int64)

returns a `Vector` of `Int64`s giving the class index of each atom
"""
function getAtomClasses(nCNA::normalCNAProfile, classes::normalCNAProfile)
    natoms = length(nCNA)
    nclasses = length(classes)

    # fill vector of size natoms with value <number-of-classes>+1 (should be 64, i.e. default class is "unclassified"i
    c = fill(nclasses, natoms)

    # determine and set class of each atom. If a class cannot be determined,
    # its default value is already "unclassified" (64 or nclasses+1).
    for i in 1:natoms
        for j in 1:nclasses
            if nCNA[i] == classes[j]
                c[i] = j
                break
            end
        end
    end
    return c
end

"""
	getFrequencyClassVector(atomClasses::Vector{Int64}, nClasses::Int64, returnType::DataType)

returns a `Vector`of type `returnType` giving the frequency of occurence for each class.
`returnType` should be UInt8 if the structure is 255 atoms or less, UInt16 if 65536 atoms or less, and so on.
`atomClasses` should be taken from `getAtomClasses()`.
"""
function getFrequencyClassVector(atomClasses::Vector{Int64}, nClasses::Int64, returnType::DataType)
	freq = Vector{returnType}(undef, nClasses)
	nAtoms = length(atomClasses)
	for i in 1:nClasses
		freq[i] = count(x->x==i, atomClasses)
	end

	return freq
end

"""
	getFractionalClassVector(atomClasses::Vector{Int64}, nClasses::Int64)

returns a `Vector`of type `Float64` giving the fractional occurence of occurence for each class.
`atomClasses` should be taken from `getAtomClasses()`.
"""
function getFractionalClassVector(atomClasses::Vector{Int64}, nClasses::Int64)
	frac = Vector{Float64}(undef, nClasses)
	nAtoms = length(atomClasses)
	for i in 1:nClasses
		frac[i] = count(x->x==i, atomClasses)/nAtoms
	end

	return frac
end

"""
	CNAToString(cna::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})

Returns a human readable string of the given CNA profile (`Dict` type)
"""
function CNAToString(cna::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})
	s = ""
	for sig in keys(cna)
		s *= "($(sig[1]), $(sig[2]), $(sig[3])): $(cna[sig]) "
	end
	return s
end

"""
	CNAToString(cna::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})

Returns a human readable string of the given CNA profile (`Vector` type)
"""
function CNAToString(cna::CNAProfile)
	s = ""
	for i in 1:length(cna)
		s *= "($(cna[i].first[1]), $(cna[i].first[2]), $(cna[i].first[3])): $(cna[i].second) "
	end
	return s
end

"""
	CNAToString(cna::Dict{Tuple{UInt8, UInt8, UInt8}, UInt16})

Returns a string describing `cna` for logging (probably into CNAlog.txt).
"""
function CNAToLogString(cna::CNAProfile)
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
(i.e. from `CNAToString(cna::CNAProfile)`)and converts it to a CNA profile (type `Vector`)
where:
	ncn,nb,nl is the CNA signagure, e.g. (5, 5, 5) = "5,5,5"
	freq is the frequency 
"""
function stringToCNA(s::String)
	
	freqCNAPair = [split(x, ':') for x in split(s, ';')]
	n = length(freqCNAPair[length(freqCNAPair)]) == 1 ? length(freqCNAPair) -1 : length(freqCNAPair)
	CNA = CNAProfile(undef, n)
	
	# For all CNA signature & frequency pairs
	for i in 1:n
		cna::Tuple{UInt8, UInt8, UInt8} = Tuple(parse(UInt8, String(x)) for x in split(freqCNAPair[i][1], ','))
		freq::UInt16 = parse(UInt16, String(freqCNAPair[i][2]))
		CNA[i] = Pair(cna, freq)
	end

	return CNA
end

"""
	stringToCNA(s::String)

Returns the CNA profile (logged as in `CNAToString(cna::CNAprofile)`) of the cluster
described by `name`.
"""
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
	elseif  name == "LJ75_Ih"
		return "5,5,5:42;4,3,3:14;4,2,2:133;3,2,2:55;3,1,1:55;3,0,0:2;2,1,1:7;2,0,0:20;"
	elseif  name == "LJ98_FCC1"
		return "4,2,2:42;4,2,1:204;3,2,2:9;3,1,1:91;3,0,0:2;2,1,1:59;2,0,0:15;1,0,0:2;"
	elseif  name == "LJ98_FCC2"
		return "5,5,5:3;4,2,2:70;4,2,1:173;3,2,2:17;3,1,1:87;3,0,0:2;2,1,1:50;2,0,0:23;1,0,0:1;"
	elseif  name == "LJ98_Dh"
		return "5,5,5:4;4,2,2:85;4,2,1:143;3,2,2:15;3,1,1:101;3,0,0:18;2,1,1:35;2,0,0:27;"
	elseif  name == "LJ98_Ih"
		return "5,5,5:27;4,2,2:159;4,2,1:62;3,2,2:60;3,1,1:87;3,0,0:7;2,1,1:14;2,0,0:21;"
	end
end