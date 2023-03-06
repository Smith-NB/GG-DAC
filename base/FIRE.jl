include("Atoms.jl")
include("LJ.jl")

using LinearAlgebra

abstract type Optimizer end


"""
    FIRE

# Arguments

- `atoms::Cluster`: The cluster to optimize
- `dt::Float64`: Some parameter for FIRE. Default values from ASE used.
- `maxstep::Float64`: Some parameter for FIRE. Default values from ASE used.
- `dtmax::Float64`: Some parameter for FIRE. Default values from ASE used.
- `NMin::Int64`: Some parameter for FIRE. Default values from ASE used.
- `finc::Float64`: Some parameter for FIRE. Default values from ASE used.
- `fdec::Float64`: Some parameter for FIRE. Default values from ASE used.
- `astart::Float64`: Some parameter for FIRE. Default values from ASE used.
- `fa::Float64`: Some parameter for FIRE. Default values from ASE used.
- `a::Float64`: Some parameter for FIRE. Default values from ASE used.
- `downhillCheck::Bool`: Some parameter for FIRE. Default values from ASE used.
"""
mutable struct FIRE <: Optimizer
	atoms::Cluster
	dt::Float64
	maxstep::Float64
	dtmax::Float64
	NMin::Int64
	finc::Float64
	fdec::Float64
	astart::Float64
	fa::Float64
	a::Float64
	downhillCheck::Bool
end

FIRE(atoms::Cluster) = FIRE(atoms, 0.1, 0.2, 1.0, 5, 1.1, 0.5, 0.1, 0.99, 0.1, false)


"""
	run(opt::FIRE, fmax::Float64)

Optimize a structure with FIRE. 

Convergence criterion: the maximum force acting upon any atom is less than fmax.
"""
function run(opt::FIRE, fmax::Float64)
	natoms = getNAtoms(opt.atoms)
	f = opt.atoms.forces
	v = zeros(natoms, 3)
	vf = dot(f, v)
	NSteps = 0
	i = 0
	while maximum(sum(f.^2, dims=2)) > fmax^2
		i += 1
		is_uphill = false
		vf = dot(f, v)

		if vf > 0 && !is_uphill
			v = (1 - opt.a) .* v .+ opt.a .* f ./ dot(f, f)^0.5 .* dot(v, v)^0.5
			if NSteps > opt.NMin
				opt.dt = minimum([opt.dt * opt.finc, opt.dtmax])
				opt.a *= opt.fa
			end
			NSteps += 1
		else
			v .*= 0.0
			opt.a = opt.astart
			opt.dt *= opt.fdec
			NSteps = 0
		end

		v .+= opt.dt .* f
		dr = opt.dt .* v
		normdr = dot(dr, dr)^0.5

		if normdr > opt.maxstep
			dr = opt.maxstep .* dr ./ normdr
		end
		
		opt.atoms.positions .+= dr
		calculate(opt.atoms, opt.atoms.calculator)
		f = opt.atoms.forces
	end
end


