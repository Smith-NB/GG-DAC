using Documenter, DAC
push!(LOAD_PATH,"../src/")

makedocs(sitename="Divide and Conquer Algorithm",
	format=Documenter.HTML(ansicolor=true), # Enable ANSI color output
		pages = [
			"Overview" => "index.md",
			"Manual" => [
				"cluster.md",
				"clustervector.md",
				"cna.md",
				"calculators.md",
				"optimiser.md",
				"metc.md",
				"reseed.md",
				"perturber.md",
				"BHA" => "bha.md"
				]
			]
		)

