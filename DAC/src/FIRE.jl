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

FIRE() = FIRE(0.1, 0.2, 1.0, 5, 1.1, 0.5, 0.1, 0.99, 0.1, false)


"""
	run(opt::FIRE, fmax::Float64)

reset the variable parameters of FIRE to their default values. Possibly redundant,
but done for safety.
"""
function reset_params!(opt::FIRE)
	opt.dt = 0.1
	opt.a = 0.1
	return nothing
end

"""
	run(opt::FIRE, fmax::Float64)

Optimize a structure with FIRE. 

Convergence criterion: the maximum force acting upon any atom is less than fmax.
"""
function optimize!(opt::FIRE, atoms::Atoms, fmax::Float64)
	natoms = getNAtoms(atoms)
	#calculateForces!(atoms, atoms.calculator)
	f = getForces!(atoms)
	v = zeros(3, natoms)
	#vf = dot(f, v)
	NSteps = 0
	i = 0
	while maximum(sum(f.^2, dims=1)) > fmax^2
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
		
		atoms.positions .+= dr'
		#calculateForces!(atoms, atoms.calculator)
		f = getForces!(atoms)
		#calculateEnergy!(atoms, atoms.calculator)
		#rintln(i, " ", atoms.energy, " ", maximum(sum(f.^2, dims=1)).^0.5)
		#f = atoms.forces
	end
	#reset_params!(opt)
	atoms.validCNA = false
	atoms.validEnergies = false
	atoms.validStresses = false
	atoms.validForces = true
	return i
end


