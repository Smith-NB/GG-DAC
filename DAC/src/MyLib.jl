"""
	removefirst!(a::Array, item::Any)

Takes an array and attempts to remove the first occurance of
the given item from it.
"""
function removefirst!(a::Array, item::Any)
	deleteat!(a, findfirst(x->x==item, a))
end

"""
	removeall!(a::Array, item::Any)

Takes an array and attempts to remove the all occurances of
the given item from it.
"""
function removeall!(a::Array, item::Any)
	deleteat!(a, findall(x->x==item, a))
end

function printsep(sep, x)
    print(x, sep)
end

function printsep(sep, xs...)
    for x in xs
        printsep(sep, x)
    end
    return nothing
end

function printsep(sep, x)
    print(x, sep)
end

function printlnsep(sep, xs...)
    for x in xs
        printsep(sep, x)
    end
    print("\n")
    return nothing
end

function binarySearch(A::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}, n::Int, T::Tuple{Int, Int, Int})
    L = 1
    R = n
    while L <= R
        m = trunc(Int, L + (R - L) / 2)
        if A[m].first > T
            L = m+1
        elseif A[m].first < T
            R = m -1
        else
            return m
        end
    end

    return -1
end

function binarySearch(A::Vector{Pair{Tuple{UInt8, UInt8, UInt8}, UInt16}}, n::Int, T::Tuple{UInt8, UInt8, UInt8})
    L = 1
    R = n
    while L <= R
        m = trunc(Int, L + (R - L) / 2)
        if A[m].first > T
            L = m+1
        elseif A[m].first < T
            R = m -1
        else
            return m
        end
    end

    return -1
end

"""
    sphericalToCartesian(ρ::Float64, θ::Float64, φ::Float64)

Takes a radius `ρ` and two spherical coordinates, `φ` and `θ` and converts them
into 3D Cartesian coordinates, returning a Vector of these coordinates.
"""
function sphericalToCartesian(ρ::Float64, θ::Float64, φ::Float64)
    φ *= π/180
    θ *= π/180

    x = ρ * sin(φ) * cos(θ)
    y = ρ * sin(φ) * sin(θ)
    z = ρ * cos(φ)

    return [x, y, z]
end