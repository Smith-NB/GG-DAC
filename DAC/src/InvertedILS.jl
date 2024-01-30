using BenchmarkTools
using PyPlot
import Random
using MultivariateStats
Random.seed!(1)

f(i::Int64, nCols::Int64) = ((i-1)÷nCols+1, (i-1)%nCols+1)
F(a::Int64, b::Int64, nCols::Int64) = (a-1)*nCols + b

"""
    inv_getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, nDims::Int64)

Calculate the Cartesian distances in n-dim space between point `X` and the points `Y`, each pointing to values in `data`. Store the resulting 
distance matrix in `D`. `lab` and `unl` each point to the positions of `X` and `Y` respectively within `D`, which are different to in `data`.   
"""
function inv_getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, nDims::Int64)
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
    inv_getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, lab::Int64, unl::Vector{Int64}, X::Int64, Y::Vector{Int64}, closestUnl::Vector{Int64}, closestUnlDist::Vector{Float64}, nDims::Int64)

Calculate the Cartesian distances in n-dim space between point `X` and the points `Y`, each pointing to values in `data`. Store the resulting 
distance matrix in `D`. `lab` and `unl` each point to the positions of `X` and `Y` respectively within `D`, which are different to in `data`.   
"""
function inv_getDistances!(data::Matrix{Float64}, D::Matrix{Float64}, nLab::Int64, unlPerm::Vector{Int64}, X::Int64, Y::Vector{Int64}, closestLab::Vector{Int64}, closestLabDist::Vector{Float64},  nDims::Int64)
    for j in 1:length(Y)
        r::Float64 = 0.0
        for k in 1:nDims
            r += (data[k, X] - data[k, Y[j]])^2
        end
        
        if r < closestLabDist[unlPerm[j]]
            closestLabDist[unlPerm[j]] = r
            closestLab[unlPerm[j]] = nLab
        end
        D[unlPerm[j], nLab] = r
    end
    
    return D
end

"""
    incrementIfGorE(a::Int32, b::Int64) 

Returns the decrement of `a` if `a` > `b`. Intended for use in `broadcast!` function
"""
inv_decrementIfGThan(a::Int64, b::Int64) = a > b ? a-1 : a


function inv_getMaxFromClosestLab!(D::Matrix{Float64}, closestLab::Vector{Int64}, closestLabDist::Vector{Float64}, unlPerm::Vector{Int64}, toUpdate::Vector{Int64}, nUnl::Int64)
   #v::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true} = Base.view(closestLabDist, unlPerm)
   v = Base.view(closestLabDist, unlPerm)
   max::Float64, unlIndex::Int64 = findmax(v) # min is the smallest distance. labIndex gives the labelled point involved.
   labIndex = closestLab[unlPerm[unlIndex]]
   #println("closestLab $(Base.view(closestLab, unlPerm))")
   #println("closestLabDist $v")
   #println("unlPerm$unlPerm")
   #:println("max unlindex labindex $max $unlIndex $labIndex")
   nToUpdate::Int64 = 0
   
   for i in 1:length(v)
       if closestLab[unlPerm[i]] == labIndex
           nToUpdate += 1
           toUpdate[nToUpdate] = unlPerm[i]
       end
   end
    
   broadcast!(inv_decrementIfGThan, closestLab, closestLab, unlIndex)
   return max, (unlIndex, labIndex), nToUpdate
end

function invILS(data::Matrix{Float64}, labels::Vector{Int64}, iterative::Bool)
    
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
    closestUnlDist = Vector{Float64}(undef, nSamples-1)            # Values for distances of above.

    closestLab = Vector{Int64}(undef, nSamples-1)
    closestLabDist = fill(Inf, nSamples-1)

    toUpdate = Vector{Int64}(undef, nSamples-1)

    nCols = size(data)[2]-1
    nLab = length(labelled)
    nUnl = length(unlabelled)

    while length(unlabelled) > 0
        
        # get distances between each labelled point to each unlabelled point
        
        inv_getDistances!(data, D, nLab, unlPerm, labelled[nLab], unlabelled, closestLab, closestLabDist, nDims)        
        Ri[i], index, nToUpdate = inv_getMaxFromClosestLab!(D, closestLab, closestLabDist, unlPerm, toUpdate, nUnl)
        
        for j in 1:nLab
            #D[index[1], j] = Inf
        end

        #println("lab $(length(labelled)) unl $(length(unlabelled)) index $index")
        labels[unlabelled[index[1]]] = labelled[index[2]]
        iterationLabelledAt[unlabelled[index[1]]-1] = i
        
        push!(labelled, popat!(unlabelled, index[1]))       # add newly labelled point to list, and remove from unlabelled
        popat!(unlPerm, index[1])                           # remove the pointer to D of the unlabelled datapoint
        
        # toUpdate gives all indices in closestUnl/closestUnlDist that need updating        
        if length(unlabelled) != 0
            for j in 1:nToUpdate
                #v::SubArray{Float64, 1, Matrix{Float64}, Tuple{Vector{Int64}, Int64}, false} = Base.view(D, toUpdate[j], 1:nLab)
                v = Base.view(D, toUpdate[j], 1:nLab)
                dist = Inf
                loc = 0
                for k in 1:length(v)
                    #println("toUpdate $(toUpdate[j]) v[k] $(v[k])")
                    if v[k] < dist
                        dist = v[k]
                        loc = k
                    end
                end
                
                closestLabDist[toUpdate[j]] = dist
                closestLab[toUpdate[j]] = loc
            end
        end
        
        i += 1
        nLab += 1
        nUnl -= 1
    end
    return Ri, iterationLabelledAt
end


