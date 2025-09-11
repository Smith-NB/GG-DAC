# Displacement Operators

```@meta
CurrentModule = DAC
```

The following displacement operators are available for perturbing a `Cluster`.

```@docs
cartesianDisplacement(coords::Matrix{Float64}, dr::Float64)
geometricCentreDisplacement(coords::Matrix{Float64}, alphaMin::Float64, alphaMax::Float64, w::Float64)
perturbClusterSurface(coords::Matrix{Float64}, nAtomsToMove::Int64, rCut::Float64)
```

# Building Composite Operators

The [`BasinHopper`](@ref) contains a `perturber` field, this being a function called to perturb the structure of the `Cluster` for each attempted hop. Even if only a single displacement operator is being used, one must define a custom function for this. Below is an example of such a function which uses two different operators, selected at random. The function must return the new coordinates, and a `String` describing which operator was used.

```julia
dr = 1.0
nAtomsToMove = 1
rcut = 3.4765
function perturber(coords::Matrix{Float64})
    r = rand()
    if r < 0.75 # this is the likelihood of the CDO being used for perturbation
        return DC.cartesianDisplacement(coords, dr), "CDO perturbation"
    end
    
    return DC.perturbClusterSurface(coords, nAtomsToMove, rcut), "SAO perturbation" 
end
```