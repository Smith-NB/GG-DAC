# Cluster

```@meta
CurrentModule = DAC
```

# The `Cluster` type

`Cluster` is a `struct` (similar to an `Object` in Python) that stores information about a nanoparticle. 

```@docs
Cluster
```
\

# Creating `Cluster` instnaces

\
Typically, it is impractical to create an instance of `Cluster` directly, so instead there are various methods capable of creating it. These are described below.
\
\
 
```@docs
Cluster(formula::Dict{String, Int64})
Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64})
Cluster(formula::Dict{String, Int64}, positions::Matrix{Float64}, cell::Matrix{Float64})
```
\

# Reading .xyz files

```@docs
read_xyz(filename::String)
read_xyzs(filename::String, formula::Dict{String, Int64})
```
\

# Writing .xyz files

```@docs
write_xyz(filename::String, atoms::Cluster)
write_xyz(filename::String, positions::Matrix{Float64}, formula::Dict{String, Int64}, cell::Float64)
write_xyz(filename::String, atoms::Cluster, tags::Vector{Int64})
write_xyz(filename::String, positions::Matrix{Float64}, formula::Dict{String, Int64}, cell::Float64, tags::Vector{Int64})
```

# Generating a random `Cluster`

```@docs
	generateRandomSeed(formula::Dict{String, Int64}, boxLength::Number, vacuumAdd::Number, returnCoordsOnly::Bool=false)
```
