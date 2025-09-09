# CNA Profiles and the SCM

```@meta
CurrentModule = DAC
```

The DAC package can calculate total and normal CNA profiles of nanoparticles. The following aliases are used as types to store these profiles.

```@docs
CNAProfile
normalCNAProfile
```

\

# Calculating a CNA profile

Note that all CNA calculation methods take in a coordinates argument (`Matrix{Float64}`) but also have wrapper functions that allow `Cluster` types to be used as arguments.

The following methods can be used to calculate CNA profiles.

```@docs
getCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
getNormalCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
getTotalAndNormalCNAProfile(coordinates::Matrix{Float64}, rcut::Float64)
```

\

# Calculating SCM similarity.

```@docs
getCNASimilarity(x::CNAProfile, y::CNAProfile)
```