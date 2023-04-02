#include("Atoms.jl")

"""
   LJ

   # Arguments

- `epsilon::Float64`: epsilon parameter of the Lennard-Jones potential.
- `sigma::Float64`: sigma parameter of the Lennard-Jones potential.
- `rc::Float64`: Cut-off distance.
"""
struct LJ_ASAP3 <: Calculator
	_py::PyObject
end

struct LJ_ASAP3_BadAtomsException <: Exception
	type1::Type
	type2::Type
end

showerror(io::IO, e::LJ_ASAP3_BadAtomsException) = print(io, "LJ_ASAP3_BadAtomsException !: $type2 cannot calculate on $type1")

calculateEnergy!(atoms::Workhorse, calc::LJ_ASAP3) = nothing
calculateEnergies!(atoms::Workhorse, calc::LJ_ASAP3) = nothing
calculateForces!(atoms::Workhorse, calc::LJ_ASAP3) = nothing
calculateStresses!(atoms::Workhorse, calc::LJ_ASAP3) = nothing

calculateEnergy!(atoms::Cluster, calc::LJ_ASAP3) = throw(LJ_ASAP3_BadAtomsException(typeof(atoms), typeof(calc)))
calculateEnergies!(atoms::Cluster, calc::LJ_ASAP3) = throw(LJ_ASAP3_BadAtomsException(typeof(atoms), typeof(calc)))
calculateForces!(atoms::Cluster, calc::LJ_ASAP3) = throw(LJ_ASAP3_BadAtomsException(typeof(atoms), typeof(calc)))
calculateStresses!(atoms::Cluster, calc::LJ_ASAP3) = throw(LJ_ASAP3_BadAtomsException(typeof(atoms), typeof(calc)))