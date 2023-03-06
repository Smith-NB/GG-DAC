include("Atoms.jl")

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
end



"""
	calculate(atoms::Cluster, calc::LJ)

Calculates the atomic energies, forces, and stresses of a Cluster type.
"""
function calculate(atoms::Cluster, calc::LJ)
	
	natoms = getNAtoms(atoms)

	e0 = 4 * calc.epsilon * ((calc.sigma / calc.rc)^12 - (calc.sigma / calc.rc)^6)
	
	energies = zeros(Float64, natoms)
	forces = zeros(Float64, natoms, 3)
	stresses = zeros(Float64, natoms, 3, 3)

	
	for ii in 1:natoms
		distance_vectors = zeros(Float64, natoms-ii, 3)
		r2 = zeros(Float64, natoms-ii)
		for jj in ii+1:natoms
			distance_vectors[jj-ii, 1] = atoms.positions[jj, 1] - atoms.positions[ii, 1]
			distance_vectors[jj-ii, 2] = atoms.positions[jj, 2] - atoms.positions[ii, 2]
			distance_vectors[jj-ii, 3] = atoms.positions[jj, 3] - atoms.positions[ii, 3]

			r2[jj-ii] = distance_vectors[jj-ii, 1]^2 + distance_vectors[jj-ii, 2]^2 + distance_vectors[jj-ii, 3]^2
		end
		#println("dv    ", distance_vectors)
		#println("dv[:] ", distance_vectors[:])
		c6 = r2.\calc.sigma^2
		c6 = c6.^3
		c12 = c6.^2


		pairwise_energies =  4 * calc.epsilon .* (c12 .- c6) .- e0
		energies[ii] += 0.5 * sum(pairwise_energies)

		pairwise_forces = reshape(-24 * calc.epsilon .* (2 .* c12 .- c6) ./ r2, :, 1) .* distance_vectors
		forces[ii, :] .+= sum(pairwise_forces, dims=1)[:]
		stresses[ii, :, :] .+= 0.5 .* pairwise_forces' * distance_vectors

		for jj in ii+1:natoms
			energies[jj] += 0.5 * pairwise_energies[jj-ii]
			forces[jj, :] .+= -pairwise_forces[jj-ii, :]
			stresses[jj, :, :] .+= 0.5 .* pairwise_forces[jj-ii, :] * distance_vectors[jj-ii, :]'
		end

		
	end
	atoms.energy = sum(energies)
	atoms.energies = energies
	atoms.forces = forces
	atoms.stresses = stresses
end

