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
	v::Matrix{Float64}
	f2::Matrix{Float64}
	f2sum::Vector{Float64}
	dr::Matrix{Float64}
end

FIRE() = FIRE(0.1, 0.2, 0.5, 5, 1.1, 0.5, 0.1, 0.99, 0.1, false, zeros(3, 1), zeros(3, 1), zeros(1), zeros(3, 1))


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

function getForces2sum!(opt::FIRE, f::Matrix{Float64}, N::Int64)
	for i in 1:length(f)
		opt.f2[i] = f[i]*f[i]
	end
	for i in 1:N
		opt.f2sum[i] = opt.f2[1, i] + opt.f2[2, i] + opt.f2[3, i]
	end

end

function myDot(a::Matrix{Float64}, b::Matrix{Float64})
	ab::Float64 = 0.0
	for i in 1:length(a)
		ab += a[i] * b[i]
	end
	return ab
end

"""
	run(opt::FIRE, fmax::Float64)

Optimize a structure with FIRE. 

Convergence criterion: the maximum force acting upon any atom is less than fmax.
"""
function optimize!(opt::FIRE, atoms::Atoms, fmax::Float64)
	natoms = getNAtoms(atoms)
	f = getForces!(atoms)
	if size(opt.v) != (3, natoms) 
		opt.v = zeros(3, natoms)
		opt.f2 = zeros(3, natoms)
		opt.f2sum = zeros(natoms)
		opt.dr = zeros(3, natoms)
	end
	opt.dt = 0.1
	opt.a = 0.1
	getForces2sum!(opt, f, natoms)
	
	fill!(opt.v, 0.0)
	fill!(opt.dr, 0.0)

	NSteps::Int64 = 0
	n = 0
	fmax2::Float64 = fmax * fmax
	while maximum(opt.f2sum) > fmax2
		n += 1

		is_uphill::Bool = false

		vf::Float64 = myDot(f, opt.v)
		if vf > 0 && !is_uphill# && n > 2
			dotf::Float64 = sqrt(myDot(f, f))
			dotv::Float64 = sqrt(myDot(opt.v, opt.v))
			for i in 1:length(opt.v)
				opt.v[i] = (1 - opt.a) * opt.v[i] + opt.a * f[i] / dotf * dotv
			end
			
			if NSteps > opt.NMin
				opt.dt = opt.dt * opt.finc
				if opt.dtmax < opt.dt
					opt.dt = opt.dtmax
				end
				
				opt.a *= opt.fa
			end
			NSteps += 1
		else
			opt.v .*= 0.0
			opt.a = opt.astart
			opt.dt *= opt.fdec
			NSteps = 0
		end

		#opt.v .+= opt.dt .* f
		for i in 1:length(opt.v)
			opt.v[i] += opt.dt * f[i]
		end
		
		for i in 1:length(opt.v)
			opt.dr[i] = opt.dt * opt.v[i]
		end
		
		normdr::Float64 = sqrt(myDot(opt.dr, opt.dr))

		if normdr > opt.maxstep
			for i in 1:length(opt.dr)
				opt.dr[i] = opt.dr[i] * opt.maxstep / normdr
			end
			
		end

		for i in 1:length(atoms.positions)
			atoms.positions[i] += opt.dr'[i]
		end
		

		f = getForces!(atoms)
		
		getForces2sum!(opt, f, natoms)
	end
	#reset_params!(opt)
	setDistances!(atoms)
	calculateEnergy!(atoms, atoms.calculator)
	atoms.validCNA = false
	atoms.validEnergies = false
	atoms.validStresses = false
	atoms.validForces = true

	return n+1
end

function optimizeOldAndCrap!(opt::FIRE, atoms::Atoms, fmax::Float64)
	natoms = getNAtoms(atoms)
	#calculateForces!(atoms, atoms.calculator)
	f = getForces!(atoms)
	f2 = f.^2
	f2sum = sum(f2, dims=1)
	if size(opt.v) != (3, natoms)
		opt.v = zeros(3, natoms)
	end
	dr::Matrix{Float64} = zeros(3, natoms)

	fill!(opt.v, 0.0)
	#vf = dot(f, v)
	NSteps::Int64 = 0
	i = 0

	#while maximum(f2sum) > fmax^2
	while maximum(sum(f.^2, dims=1)) > fmax^2
		i += 1
		is_uphill::Bool = false
		vf::Float64 = dot(f, opt.v)
		if vf > 0 && !is_uphill
			#=
			dotf::Float64 = dot(f, f)^0.5
			dotv::Float64 = dot(opt.v, opt.v)^0.5
			for i in 1:length(opt.v)
				opt.v[i] = (1 - opt.a) * opt.v[i] + opt.a * f[i] / dotf * dotv
			end
			=#
			opt.v = (1 - opt.a) .* opt.v .+ opt.a .* f ./ dot(f, f)^0.5 .* dot(opt.v, opt.v)^0.5
			if NSteps > opt.NMin
				opt.dt = minimum([opt.dt * opt.finc, opt.dtmax])
				#=
				opt.dt = opt.dt * opt.finc
				if opt.dtmax < opt.dt
					opt.dt = opt.dtmax
				end
				=#
				opt.a *= opt.fa
			end
			NSteps += 1
		else
			opt.v .*= 0.0
			opt.a = opt.astart
			opt.dt *= opt.fdec
			NSteps = 0
		end

		opt.v .+= opt.dt .* f
		dr = opt.dt .* opt.v
		#=
		for i in 1:length(opt.v)
			dr[i] = opt.dt * opt.v[i]
		end
		=#
		normdr::Float64 = dot(dr, dr)^0.5

		if normdr > opt.maxstep
			dr = opt.maxstep .* dr ./ normdr
			#=
			for i in 1:length(dr)
				dr[i] = dr[i] * opt.maxstep / normdr
			end
			=#
		end

		#=
		for i in 1:length(atoms.positions)
			atoms.positions[i] += dr'[i]
		end
		=#
		setPositions!(atoms, getPositions(atoms) + dr')

		#atoms.positions .+= dr
		#calculateForces!(atoms, atoms.calculator)
		f = getForces!(atoms)
		#=
		for i in 1:length(f)
			f2[i] = f[i]^2
		end
		for i in 1:natoms
			f2sum[i] = f2[1, i] + f2[2, i] + f2[3, i]
		end
		=#
		#calculateEnergy!(atoms, atoms.calculator)
		#rintln(i, " ", atoms.energy, " ", maximum(sum(f.^2, dims=1)).^0.5)
		#f = atoms.forces
	end
	#reset_params!(opt)
	atoms.validCNA = false
	atoms.validEnergies = false
	atoms.validStresses = false
	atoms.validForces = true
	return i+1
end


