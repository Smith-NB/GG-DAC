# Optimiser

```@meta
CurrentModule = DAC
```

Only one force optimizer has been implemented, FIRE:

```@docs
FIRE
FIRE()
```

# Post Optimisation Tasks

One field for the [`BasinHopper`](@ref) is `postOptimisationTasks`, a function called after each optimisation to clean up the `Cluster`. There are three implemented functions, shown below with their intended use cases explained. One of these functions MUST be passed to the `BasinHopper` to use the BHA or DACA.

```@docs
standardPostOptimisationTasks!(cluster::Cluster, bh::BasinHopper)
extendedPostOptimisationTasks!(cluster::Cluster, bh::BasinHopper)
extendedPostOptimisationTasksNoPCA!(cluster::Cluster, bh::BasinHopper)
```