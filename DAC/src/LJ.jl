#include("Atoms.jl")

"""
   LJ

   # Arguments

- `epsilon::Float64`: epsilon parameter of the Lennard-Jones potential.
- `sigma::Float64`: sigma parameter of the Lennard-Jones potential.
- `rc::Float64`: Cut-off distance.
"""
struct LJ <: Calculator
	epsilon::Float64
	sigma::Float64
	rc::Float64
	distance_vectors::Matrix{Float64}
	r2::Vector{Float64}
	c6::Vector{Float64}
	c12::Vector{Float64}
end



function LJ(epsilon::Number, sigma::Number, rc::Number, N::Int64)
	LJ(epsilon, sigma, rc, zeros(Float64, 3, N), zeros(Float64, N), zeros(Float64, N), zeros(Float64, N))
end

function calculateForces!(atoms::Cluster, calc::LJ, workSpace::Matrix{Float64})
	natoms = getNAtoms(atoms)
	fill!(atoms.forces, 0)

	for n in 1:natoms-1
		#fill!(calc.distance_vectors, 0)
		#fill!(calc.r2, 0)
		for i in n+1:natoms
			calc.distance_vectors[1, i] = atoms.positions[i, 1] - atoms.positions[n, 1]
			calc.distance_vectors[2, i] = atoms.positions[i, 2] - atoms.positions[n, 2]
			calc.distance_vectors[3, i] = atoms.positions[i, 3] - atoms.positions[n, 3]

			calc.r2[i] = calc.distance_vectors[1, i]*calc.distance_vectors[1, i] + calc.distance_vectors[2, i]*calc.distance_vectors[2, i] + calc.distance_vectors[3, i]*calc.distance_vectors[3, i]

			dist4 = calc.r2[i]*calc.r2[i]
			dist8 = dist4*dist4
			dist14 = dist8*dist4*calc.r2[i]

			dV = -24 * calc.epsilon*(2*calc.sigma^12/dist14 - calc.sigma^6/dist8)

			dVDist = dV * calc.distance_vectors[1, i]
			atoms.forces[1, n] += dVDist
			atoms.forces[1, i] -= dVDist

			dVDist = dV * calc.distance_vectors[2, i]
			atoms.forces[2, n] += dVDist
			atoms.forces[2, i] -= dVDist

			dVDist = dV * calc.distance_vectors[3, i]
			atoms.forces[3, n] += dVDist
			atoms.forces[3, i] -= dVDist
		end
	end
	atoms.validForces = true
end

function calculateForces!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)
	fill!(atoms.forces, 0)
	#f = zeros(3, natoms)
	#dVDist::Float64 = 0.0
	#dV::Float64 = 0.0
	
	#  1.745 s (57281724 allocations: 882.92 MiB)
	for n in 1:natoms-1
		#fill!(calc.distance_vectors, 0)
		#fill!(calc.r2, 0)
		for i in n+1:natoms
			calc.distance_vectors[1, i] = atoms.positions[i, 1] - atoms.positions[n, 1]
			calc.distance_vectors[2, i] = atoms.positions[i, 2] - atoms.positions[n, 2]
			calc.distance_vectors[3, i] = atoms.positions[i, 3] - atoms.positions[n, 3]

			calc.r2[i] = calc.distance_vectors[1, i]*calc.distance_vectors[1, i] + calc.distance_vectors[2, i]*calc.distance_vectors[2, i] + calc.distance_vectors[3, i]*calc.distance_vectors[3, i]

			dist4 = calc.r2[i]*calc.r2[i]
			dist8 = dist4*dist4
			dist14 = dist8*dist4*calc.r2[i]

			dV::Float64 = -24 * calc.epsilon*(2*calc.sigma^12/dist14 - calc.sigma^6/dist8)

			dVDist::Float64 = dV * calc.distance_vectors[1, i]
			atoms.forces[1, n] += dVDist
			atoms.forces[1, i] -= dVDist

			dVDist = dV * calc.distance_vectors[2, i]
			atoms.forces[2, n] += dVDist
			atoms.forces[2, i] -= dVDist

			dVDist = dV * calc.distance_vectors[3, i]
			atoms.forces[3, n] += dVDist
			atoms.forces[3, i] -= dVDist

			#println()
		end
	end
	#atoms.forces = f
	atoms.validForces = true
end

function calculateForcesOld!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)

	#e0 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	
	#fill!(atoms.energies, 0)
	fill!(atoms.forces, 0)
	#fill!(atoms.stresses, 0)
	for ii in 1:natoms-1
		fill!(calc.distance_vectors, 0)
		fill!(calc.r2, 0)
		for jj in ii+1:natoms
			calc.distance_vectors[1, jj-ii] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			calc.distance_vectors[2, jj-ii] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			calc.distance_vectors[3, jj-ii] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			calc.r2[jj-ii] = calc.distance_vectors[1, jj-ii]^2 + calc.distance_vectors[2, jj-ii]^2 + calc.distance_vectors[3, jj-ii]^2
		end

		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		c12 = c6.^2
		
		pairwise_forces = reshape(-24 * calc.epsilon .* (2 .* c12 .- c6) ./ calc.r2[1:natoms-ii], 1, :) .* calc.distance_vectors[:, 1:natoms-ii]
		atoms.forces[:, ii] .+= sum(pairwise_forces, dims=2)[:]

		for jj in ii+1:natoms
			atoms.forces[:, jj] .+= -pairwise_forces[:, jj-ii]
		end

	end
	atoms.validForces = true
end

function calculateEnergy!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)

	e0::Float64 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	fill!(atoms.energies, 0.0)
	#fill!(calc.c6, 0.0)
	#fill!(calc.c12, 0.0)
	#fill!(atoms.forces, 0)
	#fill!(atoms.stresses, 0)
	for ii in 1:natoms-1
		#fill!(calc.distance_vectors, 0)
		#fill!(calc.r2, 0)
		for jj in ii+1:natoms
			#=
			calc.distance_vectors[1, jj-ii] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			calc.distance_vectors[2, jj-ii] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			calc.distance_vectors[3, jj-ii] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			calc.r2[jj-ii] = 			calc.r2[jj-ii] = calc.distance_vectors[1, jj-ii]^2 + calc.distance_vectors[2, jj-ii]^2 + calc.distance_vectors[3, jj-ii]^2

			calc.c6[jj-ii] = (calc.sigma^2 / calc.r2[jj-ii])^3
			calc.c12[jj-ii] = calc.c6[jj-ii]^2
			pairwise_energy = (4 * calc.epsilon * (calc.c12[jj-ii] - calc.c6[jj-ii]) - e0) * 0.5
			=#
			dv1::Float64 = atoms.positions[jj, 1] - atoms.positions[ii, 1]
 			dv2::Float64 = atoms.positions[jj, 2] - atoms.positions[ii, 2]
 			dv3::Float64 = atoms.positions[jj, 3] - atoms.positions[ii, 3]
 			r2::Float64 = dv1^2 + dv2^2 + dv3^2

 			c6::Float64 = (calc.sigma^2 / r2)^3			

			pairwise_energy::Float64 = (4 * calc.epsilon * (c6^2 - c6) - e0) * 0.5
			atoms.energies[ii] += pairwise_energy
			atoms.energies[jj] += pairwise_energy
		end


		#=
		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		println(size(c6), " $ii")
		c12 = c6.^2


		pairwise_energies =  4 * calc.epsilon .* (c12 .- c6) .- e0
		atoms.energies[ii] += 0.5 * sum(pairwise_energies)

		for jj in ii+1:natoms
			atoms.energies[jj] += 0.5 * pairwise_energies[jj-ii]
		end
		=#

	end
	atoms.energy = sum(atoms.energies)
	atoms.validEnergies = true
end

function calculateEnergyOld!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)

	e0::Float64 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	
	fill!(atoms.energies, 0)
	#fill!(atoms.forces, 0)
	#fill!(atoms.stresses, 0)
	for ii in 1:natoms-1
		fill!(calc.distance_vectors, 0)
		fill!(calc.r2, 0)
		for jj in ii+1:natoms
			calc.distance_vectors[1, jj-ii] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			calc.distance_vectors[2, jj-ii] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			calc.distance_vectors[3, jj-ii] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			calc.r2[jj-ii] = calc.distance_vectors[1, jj-ii]^2 + calc.distance_vectors[2, jj-ii]^2 + calc.distance_vectors[3, jj-ii]^2
		end


		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		println(c6)
		c12 = c6.^2


		pairwise_energies =  4 * calc.epsilon .* (c12 .- c6) .- e0
		atoms.energies[ii] += 0.5 * sum(pairwise_energies)

		for jj in ii+1:natoms
			atoms.energies[jj] += 0.5 * pairwise_energies[jj-ii]
		end

	end
	atoms.energy = sum(atoms.energies)
	atoms.validEnergies = true
end

function calculateEnergyAndForces!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)

	e0 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	
	fill!(atoms.energies, 0)
	fill!(atoms.forces, 0)
	fill!(atoms.stresses, 0)
	for ii in 1:natoms-1
		fill!(calc.distance_vectors, 0)
		fill!(calc.r2, 0)
		#distance_vectors = zeros(Float64, 3, natoms-ii)
		#r2 = zeros(Float64, natoms-ii)
		for jj in ii+1:natoms
			calc.distance_vectors[1, jj-ii] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			calc.distance_vectors[2, jj-ii] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			calc.distance_vectors[3, jj-ii] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			calc.r2[jj-ii] = calc.distance_vectors[1, jj-ii]^2 + calc.distance_vectors[2, jj-ii]^2 + calc.distance_vectors[3, jj-ii]^2
		end


		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		c12 = c6.^2


		pairwise_energies =  4 * calc.epsilon .* (c12 .- c6) .- e0
		atoms.energies[ii] += 0.5 * sum(pairwise_energies)
			
		pairwise_forces = reshape(-24 * calc.epsilon .* (2 .* c12 .- c6) ./ calc.r2[1:natoms-ii], 1, :) .* calc.distance_vectors[:, 1:natoms-ii]
		atoms.forces[:, ii] .+= sum(pairwise_forces, dims=2)[:]
		#atoms.stresses[:, :, ii] .+= 0.5 .* pairwise_forces * calc.distance_vectors[:, 1:natoms-ii]'

		for jj in ii+1:natoms
			atoms.energies[jj] += 0.5 * pairwise_energies[jj-ii]
			atoms.forces[:, jj] .+= -pairwise_forces[:, jj-ii]
			#atoms.stresses[:, :, jj] .+= 0.5 .* pairwise_forces[:, jj-ii] * calc.distance_vectors[:, jj-ii]'
		end

	end
	atoms.energy = sum(atoms.energies)
	#atoms.energies = energies
	#atoms.forces = forces
	#atoms.stresses = stresses
	atoms.validEnergies = true
	atoms.validForces = true

end

"""
	calculate(atoms::Cluster, calc::LJ)

Calculates the atomistic energies, forces, and stresses of a Cluster type, using
the Lennard Jones potential.
"""
function calculate!(atoms::Cluster, calc::LJ)
	natoms = getNAtoms(atoms)

	e0 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	
	#atoms.energies = zeros(Float64, natoms)
	#atoms.forces = zeros(Float64, 3, natoms)
	#atoms.stresses = zeros(Float64, 3, 3, natoms)
	fill!(atoms.energies, 0)
	fill!(atoms.forces, 0)
	fill!(atoms.stresses, 0)
	for ii in 1:natoms-1
		fill!(calc.distance_vectors, 0)
		fill!(calc.r2, 0)
		#distance_vectors = zeros(Float64, 3, natoms-ii)
		#r2 = zeros(Float64, natoms-ii)
		for jj in ii+1:natoms
			calc.distance_vectors[1, jj-ii] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			calc.distance_vectors[2, jj-ii] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			calc.distance_vectors[3, jj-ii] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			calc.r2[jj-ii] = calc.distance_vectors[1, jj-ii]^2 + calc.distance_vectors[2, jj-ii]^2 + calc.distance_vectors[3, jj-ii]^2
		end


		c6 = (calc.r2[1:natoms-ii].\calc.sigma^2).^3
		c12 = c6.^2


		pairwise_energies =  4 * calc.epsilon .* (c12 .- c6) .- e0
		atoms.energies[ii] += 0.5 * sum(pairwise_energies)
			
		pairwise_forces = reshape(-24 * calc.epsilon .* (2 .* c12 .- c6) ./ calc.r2[1:natoms-ii], 1, :) .* calc.distance_vectors[:, 1:natoms-ii]
		atoms.forces[:, ii] .+= sum(pairwise_forces, dims=2)[:]
		atoms.stresses[:, :, ii] .+= 0.5 .* pairwise_forces * calc.distance_vectors[:, 1:natoms-ii]'

		for jj in ii+1:natoms
			atoms.energies[jj] += 0.5 * pairwise_energies[jj-ii]
			atoms.forces[:, jj] .+= -pairwise_forces[:, jj-ii]
			atoms.stresses[:, :, jj] .+= 0.5 .* pairwise_forces[:, jj-ii] * calc.distance_vectors[:, jj-ii]'
		end

	end
	atoms.energy = sum(atoms.energies)
	#atoms.energies = energies
	#atoms.forces = forces
	#atoms.stresses = stresses
	atoms.validEnergies = true
	atoms.validForces = true
	atoms.validStresses = true
end

calculate!(atoms::Cluster) = calculate!(atoms, atoms.calculator)
