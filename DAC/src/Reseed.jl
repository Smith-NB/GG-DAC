abstract type Reseeder end

#=============================================================================#
#===============================NewLESReseeder================================#
#=============================================================================#

mutable struct NewLESReseeder <: Reseeder
	reseedPeriod::Int64
	hopsToReseed::Int64
	reseedEnergyToBeat::Float64
	getReseedStructure::Function
	args::Vector{Any}
end

function timeToReseed!(r::NewLESReseeder)
	if r.hopsToReseed <= 0
		resetHopsToReseed!(r)
		return true
	end

	return false

end

function checkNewlyAcceptedStructure!(r::NewLESReseeder, newCluster::Cluster)
	if getEnergy(newCluster) < r.reseedEnergyToBeat
		resetHopsToReseed!(r)
		r.reseedEnergyToBeat = getEnergy(newCluster)
	end
end

updateHopsToReseed!(r::NewLESReseeder) = r.hopsToReseed -= 1

resetHopsToReseed!(r::NewLESReseeder) = r.hopsToReseed, r.reseedEnergyToBeat = r.reseedPeriod, Inf

function getReseedPeriod(r::NewLESReseeder) return r.reseedPeriod end

setHopsToReseed!(r::NewLESReseeder, hopsToReseed::Int64) = r.hopsToReseed = hopsToReseed

setReseedEnergyToBeat!(r::NewLESReseeder, energy::Float64) = r.reseedEnergyToBeat = energy

#=============================================================================#
#===============================ReseedDisabled================================#
#=============================================================================#

struct ReseedDisabled <: Reseeder end

function timeToReseed(r::ReseedDisabled) return false end

function checkNewlyAcceptedStructure(r::ReseedDisabled, newCluster::Cluster) return nothing end

function updateHopsToReseed(r::ReseedDisabled) return nothing end

function resetHopsToReseed(r::ReseedDisabled) return nothing end

function getReseedPeriod(r::ReseedDisabled) return -1 end

function setReseedPeriod!(r::ReseedDisabled, reseedPeriod::Int64) end

function setReseedEnergyToBeat!(r::ReseedDisabled, energy::Float64) end
