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

function binarySearch(A::CNAProfile, n::Int64, T::Tuple{Int, Int, Int})
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

    return -(L)
end

function binarySearch(A::CNAProfile, n::Int64, T::Tuple{UInt8, UInt8, UInt8})
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

    return -(L)
end

function binarySearch(A::Vector{Int64}, n::Int64, T::Int64)
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

    return -(L)
end


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

    return -(L)
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

function cartesianToSpherical(x::Float64, y::Float64, z::Float64)
    #x += 0.000001
    #y += 0.000001
    #z += 0.000001
    ρ = sqrt(x^2 + y^2 + z^2)
    θ = acos(z/ρ)
    φ = sign(y) * acos(x / sqrt(x^2 + y^2))
    if θ > 3.14
        θ = 0
    end
    if abs(x) < 0.01 && abs(y) < 0.01 
        φ = 0
    end
    if isnan(φ)
        φ = 0
    end
    if isnan(θ)
        θ = 0
    end
    println([ρ, θ, φ])
    println(θ == NaN)

    return [ρ, θ*180/pi, φ*180/pi]

end