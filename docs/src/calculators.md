# Calculators

```@meta
CurrentModule = DAC
```

Two `Calculator` types are implemented in DAC, described below. New `Calculator` types can be defined, and need only implement two methods: `calculateForces!` and `calculateEnergies!`. The requirements are as follows.

\

# `calculateForces!`

The method `calculateForces!(atoms::Cluster, calc::Calculator)` must set `Cluster.forces` to the force experienced by each atom in `atoms`. `Cluster.validForces` should also be set to `true`. This method should return `nothing`.

\

# `calculateEnergies!`

The method `calculateEnergies!(atoms::Cluster, calc::Calculator)` must set `Cluster.energies` to the energy experienced by each atom in `atoms`. Additionall, `atoms.energy` should be set to the final sim of `atoms.energies`. `Cluster.validEnergies` should also be set to `true`. This method should return `nothing`.

\

# Lennard-Jones

```@docs
LJ
LJ(epsilon::Float64, sigma::Float64, rc::Float64, N::Int64)
```

\

# RGL (or Gupta)

```@docs
RGL
RGL(A::Float64, p::Float64, q::Float64, r0::Float64, xi::Float64, natoms::Int64)
```