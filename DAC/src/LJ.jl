#include("Atoms.jl")

"""
   LJ

   # Arguments

- `epsilon::Float64`: epsilon parameter of the Lennard-Jones potential.
- `sigma::Float64`: sigma parameter of the Lennard-Jones potential.
- `rc::Float64`: Cut-off distance.
"""
struct LJ <: Calculator
	epsilon::Number
	sigma::Number
	rc::Number
	distance_vectors::Matrix{Float64}
	r2::Vector{Float64}
	c6::Vector{Float64}
	c12::Vector{Float64}
end


function LJ(epsilon::Number, sigma::Number, rc::Number, N::Int64)
	LJ(epsilon, sigma, rc, zeros(Float64, 3, N), zeros(Float64, N-1), zeros(Float64, N-1), zeros(Float64, N-1))
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

	#calc.distance_vectors = zeros(Float64, 3, natoms)
	#calc.r2 = zeros(Float64, natoms)
	for ii in 1:natoms-1
		#fill!(calc.distance_vectors, 0)
		#fill!(calc.r2, 0)
		#distance_vectors = zeros(Float64, natoms-ii, 3)
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
end

calculate!(atoms::Cluster) = calculate!(atoms, atoms.calculator)
