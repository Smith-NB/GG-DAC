# ClusterVector

```@meta
CurrentModule = DAC
```

When the BHA or DACA finish running, they will store all unique nanoparticle structures discovered to a file called `clusterVector.jld2`. To access these structures again, one will need to create a Julia script to read this file. The command below will load this back into Julia runtime.


```julia-repl
julia> clusterVector = DAC.jldopen("clusterVector.jld2")["clusterVector"]
```



# The `ClusterCompressed` type

```@docs
ClusterCompressed
```
\

# `ClusterVector` and `ClusterVectorWithML`

```@docs
ClusterVector
ClusterVectorWithML
```