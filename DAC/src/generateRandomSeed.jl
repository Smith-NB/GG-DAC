function isClusterCoherent(clusterCoords::Matrix{Float64}, maxDistance::Number)

	N = trunc(Int, length(clusterCoords)/3)
	neighboursList = getNeighboursList(clusterCoords, maxDistance)
	atomsNotReached = [i for i in 1:N]
	clusterPaths = [1]
	while length(clusterPaths) != 0
		atomToExplore = pop!(clusterPaths)
		index = findfirst([i==atomToExplore for i in atomsNotReached])
		if index != nothing
			deleteat!(atomsNotReached, index)
		end
		append!(clusterPaths, neighboursList[atomToExplore])
		while length(neighboursList[atomToExplore]) != 0
			atomBondedTo = neighboursList[atomToExplore][1]
			removefirst!(neighboursList[atomBondedTo], atomToExplore)
			removefirst!(neighboursList[atomToExplore], atomBondedTo)
		end
	end

	return length(atomsNotReached) == 0
end


function generateRandomSeed(formula::Dict{String, Int64}, boxLength::Number, vacuumAdd::Number, returnCoordsOnly::Bool=false)
	# get number of atoms
	N = 0
	for key in keys(formula)
		N += formula[key]
	end


	clusterCoords = zeros(Float64, N, 3)

	while true
		# Generate N atomic coordinates.
		for i in 1:N
			# Generate xyz coords for new atom and check if it is too close to other atoms.
			atomCoords = rand(Float64, 3)
			while true
				atomClash = false
				clusterCoords[i, 1:3] = rand(Float64, 3) * boxLength

				for j in 1:i-1
					if abs(clusterCoords[i, 1] - clusterCoords[j, 1]) < 0.0001 && abs(clusterCoords[i, 2] - clusterCoords[j, 2]) < 0.0001 && abs(clusterCoords[i, 3] - clusterCoords[j, 3]) < 0.0001
						atomClash = true
					end
				end

				if !atomClash
					break
				end
			end
		end

		# Check that the cluster is coherent. Restart if it isn't. 
		if isClusterCoherent(clusterCoords, 1.5)
			break
		else
			clusterCoords = zeros(Float64, N, 3)
			println("Cluster non-coherent. Trying again.")
		end
	end
	
	if returnCoordsOnly
		return clusterCoords
	else	
		cellLength = boxLength + vacuumAdd
		cell = zeros(Float64, 3, 3)
		cell[1] = cellLength; cell[5] = cellLength; cell[9] = cellLength
		return Cluster(formula, clusterCoords, cell)
	end

end

function perturbCluster(coords::Matrix{Float64}, dr::Float64)
	n = getNAtoms(coords)
	#d = [-0.05656603889502572 0.12980658943320017 0.060290983911345644; 0.3550686514429878 -0.06555336530112675 0.24440180205442533; -0.15486432997618957 0.21969807689225523 0.32139438477792526; -0.10727464246468808 -0.32967182284413665 0.25750531571476504; -0.313580204831671 -0.051923328901951483 -0.36526468482616514; -0.32537141047247176 -0.09889866529644485 -0.10987130675634821; 0.30648946597417615 0.11971258762120174 -0.2250431812127264; 0.3300962272985052 0.028945091759996447 -0.08173306381238828; 0.3223432247580048 0.10935075829141372 0.18812241338915633; -0.10080969837969588 -0.33325191272613797 -0.34078288626802034; -0.04673979517105034 -0.3713698409690649 -0.1490498508392964; 0.06498428815602501 0.0021591511038631774 0.19222540622785741; -0.28001419544715145 -0.36954934826170566 0.38296135404391696; -0.04673330002795169 -0.014489369183458934 -0.04278279670856833; 0.15974933896579754 0.24065540328142357 -0.06721624117105672; -0.17557667381476955 0.04927332258614019 -0.278046189433469; 0.26991764711203875 -0.10877293692002708 -0.3060625752242295; -0.24499332408304397 -0.24482765683406696 -0.36234549071297634; -0.2685796953746839 -0.09901396653221434 -0.32010466092617085; -0.3776662353386591 -0.3494663266636554 -0.17489042636114738; 0.014993783262024963 -0.09095545364892069 -0.3091345691686607; -0.1602151089720165 0.04586218475678941 0.08746863310283048; -0.0316637173185943 0.09919875942618016 0.14438290609972837; -0.3732770435482475 0.01843024142586085 0.11402443118913981; -0.17211398922083287 0.3523772737042931 0.15155852277638893; 0.2876285337236132 -0.32853377314132715 -0.3521592476010179; -0.21331156305506546 -0.08534985834432485 -0.16252885805506107; -0.1938876500752798 0.010281645081336156 -0.2909190440970614; 0.2708432390041535 -0.11346726904159175 0.2434789999533029; -0.33906097379982897 -0.26061220552883485 0.029056940660340126; -0.07609831639161335 -0.32093129398969333 0.3372245411204004; 0.14031058266927304 0.168297276957048 -0.11859960661168421; 0.022318827580933044 0.18928897804619105 0.28969511412436183; 0.21918758433739285 0.06180973449076737 0.06185337455702964; 0.2151543538572888 -0.2679317235997405 0.06682629676356876; 0.3069415594212062 0.12890680373547198 0.26239846000578276; -0.2516238549598055 0.07909148840140921 -0.21800954050136748; -0.3879874611019327 0.3117494999542889 -0.07045736171046686]
	#return coords + d * dr
	return coords + (rand(n, 3) * 2 .- 1) * dr
end

function perturbCluster(atoms::Cluster, dr::Float64)
	return perturbCluster(atoms.positions, dr)
end