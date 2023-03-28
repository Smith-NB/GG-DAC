#include("Atoms.jl")

"""
   LJ

   # Arguments

- `epsilon::Float64`: epsilon parameter of the Lennard-Jones potential.
- `sigma::Float64`: sigma parameter of the Lennard-Jones potential.
- `rc::Float64`: Cut-off distance.
"""
struct LJ_ASAP3 <: Calculator
	asap3Object::PyObject
end

function calculateForces!(atoms::Cluster, calc::LJ_ASAP3)

end

function calculateForces!(atoms::Workhorse, calc::LJ_ASAP3)

end