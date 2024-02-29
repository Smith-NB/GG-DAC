#include("Atoms.jl")

"""
   LJ

   # Arguments

- `epsilon::Float64`: epsilon parameter of the Lennard-Jones potential.
- `sigma::Float64`: sigma parameter of the Lennard-Jones potential.
- `rc::Float64`: Cut-off distance.
"""
struct RGL <: Calculator
	A::Float64
	p::Float64
	q::Float64
	r0::Float64
	xi::Float64
	xi2::Float64
end

function RGL(A::Float64, p::Float64, q::Float64, r0::Float64, xi::Float64)
	RGL(A, p, q, r0, xi, xi*xi)
end

function φ(k::Int64, λ::Float64, r::Matrix{Float64}, r_0::Float64, natoms::Int64)
	φ_k::Float64 = 0.0
	for l in 1:natoms
		if k == l continue end
		φ_k += exp(-λ * (r[k, l]/r_0) - 1.0)
	end

	return φ_k
end

function get_φ(atoms::Cluster, calc::RGL, natoms::Int64) 
	φ = zeros(Float64, natoms)
	for i in 1:natoms
		for j in i+1:natoms
			if i == j continue end
			r_ij::Float64 = ((atoms.positions[i, 1] - atoms.positions[j, 1])^2 + (atoms.positions[i, 2] - atoms.positions[j, 2])^2 + (atoms.positions[i, 3] - atoms.positions[j, 3])^2)^0.5

			R::Float64 = r_ij / calc.r0 - 1.0
			dq::Float64 = calc.xi2 * exp(-2.0 * calc.q * R)

			φ[i] += dq
			φ[j] += dq
		end
	end

	return φ
end

function calculateForces!(atoms::Cluster, calc::RGL)
	natoms =  getNAtoms(atoms)
	fill!(atoms.forces, 0.0)
	φ = get_φ(atoms, calc, natoms)
	for i in 1:length(φ)
		φ[i] = sqrt(φ[i])
	end

	for i in 1:natoms
		for j in i+1:natoms
			rx::Float64 = atoms.positions[i, 1] - atoms.positions[j, 1]
			ry::Float64 = atoms.positions[i, 2] - atoms.positions[j, 2]
			rz::Float64 = atoms.positions[i, 3] - atoms.positions[j, 3]
			r_ij::Float64 = sqrt(rx*rx + ry*ry + rz*rz)
			R::Float64 = r_ij / calc.r0 - 1.0

			pTerm::Float64 = -2.0 * calc.A * calc.p / calc.r0 * exp(-calc.p * R) 
			qTerm::Float64 = -calc.xi2 * calc.q / calc.r0 * exp(-2.0 * calc.q * R)
			F_ij::Float64 = ( pTerm - qTerm * ( (1.0 / φ[i]) + (1.0 / φ[j]) ) ) / r_ij

			
			Fa::Float64 = F_ij * rx
			atoms.forces[1, i] -= Fa
			atoms.forces[1, j] += Fa

			Fa = F_ij * ry
			atoms.forces[2, i] -= Fa
			atoms.forces[2, j] += Fa
			
			Fa = F_ij * rz
			atoms.forces[3, i] -= Fa
			atoms.forces[3, j] += Fa
			

		end
	end
	atoms.validForces = true

end


function calculateEnergyTime!(atoms::Cluster, calc::RGL)
	natoms = getNAtoms(atoms)
	fill!(atoms.energies, 0.0)

	sigma_p = zeros(Float64, natoms)
	sigma_q = zeros(Float64, natoms)
	t1 = 0.0
	t2 = 0.0
	t3 = 0.0
	t4 = 0.0
	t5 = 0.0
	t6 = 0.0
	t7 = 0.0
	t8 = 0.0
	t9 = 0.0
	for i in 1:natoms
		
		for j in i+1:natoms
			if i == j continue end
			t1 += @elapsed r_ij::Float64 = sqrt((atoms.positions[i, 1] - atoms.positions[j, 1])^2 + (atoms.positions[i, 2] - atoms.positions[j, 2])^2 + (atoms.positions[i, 3] - atoms.positions[j, 3])^2)

			t2 += @elapsed dx::Float64 = r_ij / calc.r0 - 1.0
			t3 += @elapsed dp::Float64 = calc.A * exp(-calc.p * dx)
			t4 += @elapsed dq::Float64 = calc.xi2 * exp(-2.0 * calc.q * dx)

			t5 += @elapsed sigma_p[i] += dp
			t6 += @elapsed sigma_q[i] += dq

			t7 += @elapsed sigma_p[j] += dp
			t8 += @elapsed sigma_q[j] += dq
		end
		t9 += @elapsed atoms.energies[i] = sigma_p[i] - sigma_q[i]^0.5
	end
	atoms.energy = sum(atoms.energies)
	atoms.validEnergies = true
	println(t1)
	println(t2)
	println(t3)
	println(t4)
	println(t5)
	println(t6)
	println(t7)
	println(t8)
	println(t9)
	return nothing

end

function calculateEnergy!(atoms::Cluster, calc::RGL)
	natoms = getNAtoms(atoms)
	fill!(atoms.energies, 0.0)

	sigma_p = zeros(Float64, natoms)
	sigma_q = zeros(Float64, natoms)
	for i in 1:natoms
		
		for j in i+1:natoms
			if i == j continue end
			r_ij::Float64 = sqrt((atoms.positions[i, 1] - atoms.positions[j, 1])^2 + (atoms.positions[i, 2] - atoms.positions[j, 2])^2 + (atoms.positions[i, 3] - atoms.positions[j, 3])^2)

			dx::Float64 = r_ij / calc.r0 - 1.0
			dp::Float64 = calc.A * exp(-calc.p * dx)
			dq::Float64 = calc.xi2 * exp(-2.0 * calc.q * dx)

			sigma_p[i] += dp
			sigma_q[i] += dq

			sigma_p[j] += dp
			sigma_q[j] += dq
		end
		atoms.energies[i] = sigma_p[i] - sqrt(sigma_q[i])
	end
	atoms.energy = sum(atoms.energies)
	atoms.validEnergies = true
	return nothing

end

