Garden Group Divide-and-Conquer Algorithm (DACA) v1.2.4

Clone the repo or otherwise download to a directory of the users choice.




Example run files are provided to run the DACA or the basin-hopping algorithm (BHA). These must be edited to point to installation directory, i.e.: `include("<install_dir>/DAC/src/DAC.jl")`
The example run files are named `BHA_example_run_file.jl` for the BHA and `DACA_example_run_fle.jl` for the DACA. These files are fully commented and should allow for the replication of all data presented in the paper.

The user must have Julia installed, version 1.8.0 is recommended.

The following packages are required to be installed in Julia as dependencies:
	
	 LinearAlgebra
	 GaussianMixtures
	 MultivariateStats
	 Random
	 HDF5
	 JLD2
	 Dates
	 BenchmarkTools

One can add the DAC package to julia manually with the following process (note the ']' key puts Julia into package mode). This will allow for one to replace the `include` line in the example run files with `using DAC` which calls on the precompiled DAC package (this is faster).
	
	> cd <install-dir>/GG-DAC
	> julia
	julia> ]
	pkg> activate DAC
	pkg> instantiate
	pkg> precompile
	

The BHA and DACA may be run from the terminal using the command: `julia <runFileName.jl>`

For the DACA, we reccomend running Julia with as many threads are there are divisions of the PES made (i.e. k). To do this, instead run: `julia -t <k> <runFileName.jl>`

The DACA requires a training dataset for the Gaussian mixture model (and fitting the Principal Component Analysis model, if enabled). If the example run file `DACA_example_run_file.jl` is used, then a file in the directory of said file, whose name begins with `explorationData` will automatically be selected to obtain this training dataset from. To generate such a file, we reccomend running the BHA (e.g. using `BHA_example_run_file.jl`) with HISTO enabled as the Metropolis Criterion, for 10,000 minimisations. This will result in a file named `clusterVector.jld2` being created, containing the structural data from the BHA run in the correct format to be used by the DACA for training. Simply copy this to the directory where the DACA is to be run, and rename the file to `explorationData.jld2`.

Two case studies are provided with example input and output in `DACA_exampleRun` and `DACA_exampleRun_noPCA`, with and without PCA enabled of runs of `Au55`, respectively.
The output files of the program are:


	`CNAlog.txt`: file containing the CNA profiles of each unique structure. Format is `<ID>=<CNAProfile>E-<Energy>`
		where `<CNAProfile` is formatted as: (`n_cn`, `n_b`, `n_l`):`freq`; ...
			`n_cn` is the number of common neighbours for the CNA signature,
			`n_b` is the number of bonds between the common neighbours,
			`n_l` is the longest chain of bonds between common neighbours.


	`DAC.out`: verbose output of the run. Each step lists the walk ID, step number (within the walk), the walker (i.e. thread) being used,
		ID of the structure encountered (negative value indicates re-encounter),
		Energy (E) of structure,
		The odds of the hop being accepted (before cluster restriction is accoutned for),
		Which cluster (division) the structure (newCluster, yes that is confusing, sorry) belong too. 
			(the above line is only listed if the hop was accepted based on energy (above probablity). it will be rejected if it belongs to the wrong cluster
		If the current step was accepted or rejected.


	`GMM_PCA.jld2`: JLD2 file storing the PCA (if enabled) and GMM models trained, see `processResults.jl` for how to access.


	`clusterVector.jld2`: JLD2 file storing unique structure positions, energies, and CNA profiles. see `processResults.jl` for how to access.

	
	`log.txt`: log file of the run. Each line lists the walkID, step number (within the walk), ID of structure (negative if re-encounter), energy, and acceptance of the step.
		From this file, MFETs are calculated. Encounter time of a given trial is taken to be the line index of `log.txt` the energy of the global minimum/target structure
		is first found on.

	
	`tasklog.txt`: used for debugging only.

In each example input/output folder, `processResults.jl` is also provided to show examples of how structural data can be retrieved from the above files, namely:
	Structure energies and similarities to a reference structure, used to construct `energy vs similarity plots`,
	the raw coordinates of structures,
	the CNA profiles of structures,
	the prevalances of different atomic classes, using `Atom-64-Class` or `Atom-80-Class`,
	the clusters (divisions) the structures belong to according to the trained Gaussian Mixture Model (and PCA if enabled),
	the classifications of the structures (i.e. FCC, TWI, ICO, PICO, AICO, DEC, or AMB).

