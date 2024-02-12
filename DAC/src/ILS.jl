using BenchmarkTools
using PyPlot
import Random
using MultivariateStats

f(i::Int64, nCols::Int64) = ((i-1)÷nCols+1, (i-1)%nCols+1)
F(a::Int64, b::Int64, nCols::Int64) = (a-1)*nCols + b

"""
    getDistances(data::Matrix{Float64}, X::Vector{Int64}, Y::Vector{Int64})

Given a given a single matrix, `data` of m datapoints in n-dimensional space, and two Vectors
`X` and `Y`, each containing a subset of these m points (and in combination containing all m points), 
calculate the distance matrix between points `X` and `Y`.
"""
function getDistances(data::Matrix{Float64}, X::Vector{Int64}, Y::Vector{Int64}, nDims::Int64)
    D = Matrix{Float64}(undef, length(X), length(Y))
    
    for i in 1:length(X)
        for j in 1:length(Y)
            r::Float64 = 0.0
            for k in 1:nDims
                r += (data[k, X[i]] - data[k, Y[j]])^2
            end
            D[i, j] = r^0.5
        end
    end
    return D
end

"""
    binarySearch(A::Vector{Float64}, n::Int64, T::Float64)

Given a sorted array of values `A`, of lenght `n`, find and return the insertion index for
the new value, `T`.
"""
function binarySearch(A::Vector{Float64}, n::Int64, T::Float64)
    L = 1
    R = n
    while L <= R
        m = trunc(Int, L + (R - L) / 2)
        if A[m] < T
            L = m+1
        elseif A[m] > T
            R = m -1
        else
            return m
        end
    end

    return L
end


"""
    getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, nDims::Int64)

Calculate the Cartesian distances in n-dim space between point `X` and the points `Y`, each pointing to values in `data`. Store the resulting 
distance matrix in `D`. `lab` and `unl` each point to the positions of `X` and `Y` respectively within `D`, which are different to in `data`.   
"""
function getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, nDims::Int64)
    for j in 1:length(Y)
        r::Float64 = 0.0
        for k in 1:nDims
            r += (data[k, X] - data[k, Y[j]])^2
        end
        D[lab, unl[j]] = r^0.5
        #D[unl[j], lab] = r^0.5
    end
 
    return D
end

"""
    getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, closestUnl::Vector{Int64}, closestDist::Vector{Float64}, nDims::Int64)

Calculate the Cartesian distances in n-dim space between point `X` and the points `Y`, each pointing to values in `data`. Store the resulting 
distance matrix in `D`. `lab` and `unl` each point to the positions of `X` and `Y` respectively within `D`, which are different to in `data`.   
"""
function getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, nLab::Int64, unlPerm::Vector{Int64}, X::Int64, Y::Vector{Int64}, closestUnl::Vector{Int64}, closestDist::Vector{Float64},  nDims::Int64)
    min::Float64, minIndex::Int64 = Inf, 0
    #min = fill(Inf, Threads.nthreads())
    #minIndex = fill(0, Threads.nthreads())
    #println("foroop")
    for j in 1:length(Y)
        r::Float64 = 0.0
        #tID = Threads.threadid()
        for k in 1:nDims
            r += (data[k, X] - data[k, Y[j]])^2
        end
        if r < min
            min = r
            minIndex = j
            #min[tID] = D[unlPerm[j], nLab]
            #minIndex[tID] = j
        end
        D[unlPerm[j], nLab] = r
    end
    
    closestUnl[nLab] = minIndex
    closestDist[nLab] = min
    
    return D
end

"""
    incrementIfGorE(a::Int32, b::Int64) 

Returns the decrement of `a` if `a` > `b`. Intended for use in `broadcast!` function
"""
decrementIfGThan(a::Int64, b::Int64) = a > b ? a-1 : a


function getMinFromClosestUnl!(D::Matrix{Float64}, closestUnl::Vector{Int64}, closestDist::Vector{Float64}, unlPerm::Vector{Int64}, toUpdate::Vector{Int64}, nLab::Int64)
   v::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} = Base.view(closestDist, 1:nLab)
   min::Float64, labIndex::Int64 = findmin(v) # min is the smallest distance. labIndex gives the labelled poiint involved.
   unlIndex = closestUnl[labIndex]
   nToUpdate::Int64 = 0
   
   for i in 1:length(closestUnl)
       if closestUnl[i] == unlIndex
           nToUpdate += 1
           toUpdate[nToUpdate] = i
       end
   end
    
   broadcast!(decrementIfGThan, closestUnl, closestUnl, unlIndex)
   return min, (unlIndex, labIndex), nToUpdate
end

function ILS(data::Matrix{Float64}, labels::Vector{Int64}, iterative::Bool)
    
    # seperate labelled and unlabelled points
    nDims, nSamples = size(data)
    
    labelled = findall(i->i!=0, labels)             # indices of all labelled points in data
    unlabelled = findall(i->i==0, labels)           # indices of all unlabelled points in data
    unlPerm = [i for i in 1:length(unlabelled)]     # Gives the row index in D for the corresponding unlabelled point in `unlabelled`.
    
    Ri = Vector{Float64}(undef, length(unlabelled))             # The distances between the closest labelled and unlabelled points at each interation
    iterationLabelledAt = Vector{Int64}(undef, length(unlabelled))   # The order points were labelled in. iterationLabelledAt[1] gives the interation the first unlabelled point was labelled at.
    
    i = 1
    D = Matrix{Float64}(undef, nSamples-1, nSamples-1)          # Distance matrix working space
    #D = zeros(nSamples-1, nSamples-1)
    closestUnl = Vector{Int64}(undef, nSamples-1)               # Gives the closest unlabelled point to each labelled point. closestUnl[1] is the closest unlabelled point to labelled[1]
    closestDist = Vector{Float64}(undef, nSamples-1)            # Values for distances of above.

    toUpdate = Vector{Int64}(undef, nSamples-1)

    nCols = size(data)[2]-1
    nLab = length(labelled)
    nUnl = length(unlabelled)

    while length(unlabelled) > 0
        
        # get distances between each labelled point to each unlabelled point
        
        getDistances!(data, D, nLab, unlPerm, labelled[nLab], unlabelled, closestUnl, closestDist, nDims)        
        Ri[i], index, nToUpdate = getMinFromClosestUnl!(D, closestUnl, closestDist, unlPerm, toUpdate, nLab)
        
        for j in 1:nLab
            D[index[1], j] = Inf
        end
        
        labels[unlabelled[index[1]]] = labelled[index[2]]
        iterationLabelledAt[unlabelled[index[1]]-1] = i
        
        push!(labelled, popat!(unlabelled, index[1]))       # add newly labelled point to list, and remove from unlabelled
        popat!(unlPerm, index[1])                           # remove the pointer to D of the unlabelled datapoint
        
        # toUpdate gives all indices in closestUnl/closestDist that need updating
        
        if length(unlabelled) != 0
            for j in 1:nToUpdate
                v::SubArray{Float64, 1, Matrix{Float64}, Tuple{Vector{Int64}, Int64}, false} = Base.view(D, unlPerm, toUpdate[j])
                dist = Inf
                loc = 0
                for k in 1:length(v)
                    if v[k] < dist
                        dist = v[k]
                        loc = k
                    end
                end
                
                closestDist[toUpdate[j]] = dist
                closestUnl[toUpdate[j]] = loc
            end
        end
        
        i += 1
        nLab += 1
        nUnl -= 1
    end
    return Ri, iterationLabelledAt
end

