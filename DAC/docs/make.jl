using Documenter, DAC
push!(LOAD_PATH,"../src/")

makedocs(sitename="Divide and Conquer Algorithm",
		pages = [
			"Overview" => "index.md",
			"Manual" => [
				"calculators.md",			
				"BHA" => "bha.md"
				]
			]
		)

