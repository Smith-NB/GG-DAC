abstract type Reseeder end

#=============================================================================#
#===============================NewLESReseeder================================#
#=============================================================================#

"""
	NewLESReseeder

Will trigger a reseed after `reseedPeriod` if no new lowest energy structure (since the last reseed) is found within that period. This period resets when a new LES is found.

# Fields

- `reseedPeriod::Int64`: The period after the latest LES discovery until a reseed is triggered.
- `hopsToReseed::Int64`: The hops until the next reseed is triggered if no LES is found.
- `reseedEnergyToBeat::Float64`: Tracks the latest LES energy.
- `getReseedStructure::Function`: Function called if a reseed is triggered to generate the new seed.
- `args::Vector{Any}`: Function arguments passed to `getReseedStructure`.

# Example

```julia
reseedPeriod = 50
reseedEnergyToBeat = Inf # any structure will have a lower energy than +.v.e infty.
formula = Dict("Au" => 55)
boxLength = 5.0
vacuumAdd = 10.0
returnCoordsOnly = true
coherencyDistance = 4.0
getReseedStructure = generateRandomSeed
args = [formula, boxLength, vacuumAdd, returnCoordsOnly, coherencyDistance]
reseeder = NewLESReseeder(reseedPeriod, reseedPeriod, reseedEnergyToBeat, 
					getReseedStructure, args)
```

"""
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

function getHopsToReseed(r::NewLESReseeder) return r.hopsToReseed end

function getReseedEnergyToBeat(r::NewLESReseeder) return r.reseedEnergyToBeat end

setHopsToReseed!(r::NewLESReseeder, hopsToReseed::Int64) = r.hopsToReseed = hopsToReseed

setReseedEnergyToBeat!(r::NewLESReseeder, energy::Float64) = r.reseedEnergyToBeat = energy


#=============================================================================#
#===============================ReseedDisabled================================#
#=============================================================================#

"""
	ReseedDisabled

Will never trigger a reseed.

# Example

```julia
reseeder = ReseedDisabled()
```
"""
struct ReseedDisabled <: Reseeder end

function timeToReseed!(r::ReseedDisabled) return false end

function checkNewlyAcceptedStructure!(r::ReseedDisabled, newCluster::Cluster) return nothing end

function updateHopsToReseed!(r::ReseedDisabled) return nothing end

function resetHopsToReseed!(r::ReseedDisabled) return nothing end

function getReseedPeriod(r::ReseedDisabled) return -1 end

function getHopsToReseed(r::ReseedDisabled) return -1 end

function getReseedEnergyToBeat(r::ReseedDisabled) return Inf end

function setHopsToReseed!(r::ReseedDisabled, reseedPeriod::Int64) end

function setReseedEnergyToBeat!(r::ReseedDisabled, energy::Float64) end

#=============================================================================#
#===============================ELimReseeder================================#
#=============================================================================#

mutable struct ELimReseeder <: Reseeder
	reseedPeriod::Int64
	hopsToReseed::Int64
	reseedEnergyToBeat::Float64
	getReseedStructure::Function
	args::Vector{Any}
	ELimBounceCounter::Int64
	eLimBounceLimit::Int64
	forceReseed::Bool
end

function timeToReseed!(r::ELimReseeder)
	if r.hopsToReseed <= 0 || r.forceReseed
		resetHopsToReseed!(r)
		r.forceReseed = false
		return true
	end

	return false

end

function checkNewlyAcceptedStructure!(r::ELimReseeder, newCluster::Cluster)
	if getEnergy(newCluster) < r.reseedEnergyToBeat
		resetHopsToReseed!(r)
		r.reseedEnergyToBeat = getEnergy(newCluster)
	end
end

updateHopsToReseed!(r::ELimReseeder) = r.hopsToReseed -= 1

resetHopsToReseed!(r::ELimReseeder) = r.hopsToReseed, r.reseedEnergyToBeat = r.reseedPeriod, Inf

function getReseedPeriod(r::ELimReseeder) return r.reseedPeriod end

function getHopsToReseed(r::ELimReseeder) return r.hopsToReseed end

function getReseedEnergyToBeat(r::ELimReseeder) return r.reseedEnergyToBeat end

setHopsToReseed!(r::ELimReseeder, hopsToReseed::Int64) = r.hopsToReseed = hopsToReseed

setReseedEnergyToBeat!(r::ELimReseeder, energy::Float64) = r.reseedEnergyToBeat = energy

function eLimCrossed!(r::ELimReseeder)
	r.ELimBounceCounter += 1 #increment counter

	# if the ELim has been bounced off too many times, set the hops to reseed to 0
	# This will force a reseed at the next call of `timeToReseed!`
	if r.ELimBounceCounter >= r.eLimBounceLimit
		r.forceReseed = true
	end
end
