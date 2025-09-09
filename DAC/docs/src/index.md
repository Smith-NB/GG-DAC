# DAC.jl Documentation

```@meta
CurrentModule = DAC
```

The DAC (Divide-and-Conquer) module implements a Divide-and-Conquer Algorithm (DACA) for atomic nanoparticle global optimisation. Two scripts are nessecary to run the DACA, one to perform exploration of the PES and gather data, and another to perform seperate, parallel walks of the subsequent divisions of the PES, made based on the exploration data.

There is currently two energy potential implemented, the Lennard-Jones (LJ) potential, and the Gupta potential (also called the RGL potential). There is currently only one energy minimiser implemented, the Fast Inertrial Relaxation Engine (FIRE).


```@docs
logCNA(io::Tuple{IO, Channel}, ID::Int64, CNA::CNAProfile, energy::Float64)
```
