
struct Workhorse <: Atoms
	_py::PyObject
end


function getFormula(atoms::Workhorse) 
	symbols = atoms._py.get_chemical_symbols()
	return Dict((i, count(==(i), symbols)) for i in unique(symbols))
end

function getNAtoms(atoms::Workhorse)
	f = getFormula(atoms)
	return sum(f[i] for i in keys(f))
end

getPositions(atoms::Workhorse) = atoms._py.get_positions()
getCell(atoms::Workhorse) = atoms._py.get_cell()
getEnergy!(atoms::Workhorse) = atoms._py.get_potential_energy()
getEnergies!(atoms::Workhorse) = atoms._py.get_potential_energies()
getForces!(atoms::Workhorse) = atoms._py.get_forces()
getStresses!(atoms::Workhorse) = atoms._py.get_stresses()
getCNA(atoms::Workhorse) = throw(:CNA_is_not_defined_for_DAC_Workhorse)
getValidCNA(atoms::Workhorse) = true
getValidEnergies(atoms::Workhorse) = true
getValidForces(atoms::Workhorse) = true
getValidStresses(atoms::Workhorse) = true
getCalculator(atoms::Workhorse) = atoms._py.get_calculator()

setPositions!(atoms::Workhorse, positions::Matrix{Float64}) = atoms._py.set_positions(positions)
setCell!(atoms::Workhorse, cell::Matrix{Float64}) = atoms._py.set_cell(cell)
setValidCNA!(atoms::Workhorse, valid::Bool) = nothing
setValidEnergies!(atoms::Workhorse, valid::Bool) = nothing
setValidForces!(atoms::Workhorse, valid::Bool) = nothing
setValidStresses!(atoms::Workhorse, valid::Bool) = nothing
setCalculator!(atoms::Workhorse, calc::_py_Calculator) = atoms._py.set_calculator(calc._py)