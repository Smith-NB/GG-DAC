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
	for i in 1:natoms
		for j in i+1:natoms
			rx::Float64 = atoms.positions[i, 1] - atoms.positions[j, 1]
			ry::Float64 = atoms.positions[i, 2] - atoms.positions[j, 2]
			rz::Float64 = atoms.positions[i, 3] - atoms.positions[j, 3]
			r_ij::Float64 = (rx^2 + ry^2 + rz^2)^0.5
			R::Float64 = r_ij / calc.r0 - 1.0

			pTerm::Float64 = -2.0 * calc.A * calc.p / calc.r0 * exp(-calc.p * R) 
			qTerm::Float64 = -calc.xi2 * calc.q / calc.r0 * exp(-2.0 * calc.q * R)
			F_ij::Float64 = ( pTerm - qTerm * ( (1.0 / φ[i]^0.5) + (1.0 / φ[j]^0.5) ) ) / r_ij

			
			atoms.forces[1, i] -= F_ij * rx
			atoms.forces[1, j] += F_ij * rx

			atoms.forces[2, i] -= F_ij * ry
			atoms.forces[2, j] += F_ij * ry

			atoms.forces[3, i] -= F_ij * rz
			atoms.forces[3, j] += F_ij * rz
			

		end
	end
	atoms.validForces = true

end


function calculateEnergy!(atoms::Cluster, calc::RGL)
	natoms = getNAtoms(atoms)
	fill!(atoms.energies, 0.0)

	sigma_p = zeros(Float64, natoms)
	sigma_q = zeros(Float64, natoms)
	for i in 1:natoms
		
		for j in i+1:natoms
			if i == j continue end
			r_ij::Float64 = ((atoms.positions[i, 1] - atoms.positions[j, 1])^2 + (atoms.positions[i, 2] - atoms.positions[j, 2])^2 + (atoms.positions[i, 3] - atoms.positions[j, 3])^2)^0.5

			dx::Float64 = r_ij / calc.r0 - 1.0
			dp::Float64 = calc.A * exp(-calc.p * dx)
			dq::Float64 = calc.xi2 * exp(-2.0 * calc.q * dx)

			sigma_p[i] += dp
			sigma_q[i] += dq

			sigma_p[j] += dp
			sigma_q[j] += dq
		end
		atoms.energies[i] = sigma_p[i] - sigma_q[i]^0.5
	end
	atoms.energy = sum(atoms.energies)
	atoms.validEnergies = true
	println("Energy sq: $sigma_q")
	return nothing

end


