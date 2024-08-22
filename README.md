Garden Group Divide-and-Conquer Algorithm (DACA) v1.2.4

Clone the repo or otherwise downloaded to a directory of the users choice.
Example run files are provided to run the DACA or the basin-hopping algorithm (BHA).
The User must have Julia installed, version 1.8.0 is recommended.

The BHA and DACA may be run using the command in the terminal `julia <runFileName.jl>`

For the DACA, we reccomend running Julia with as many threads are there are divisions of the PES made (i.e. k). To do this, instead run: `julia -t <k> <runFileName.jl>`

The DACA requires a training dataset for the Gaussian mixture model (and fitting the Principal Component Analysis model, if enabled). If the example run file `DACA_example_run_file.jl` is used, then a file in the directory of said file, whose name begins with `explorationData` will automatically be selected to obtain this training dataset from. To generate such a file, we reccomend running the BHA (e.g. using `BHA_example_run_file.jl`) with HISTO enabled as the Metropolis Criterion, for 10,000 minimisations. This will result in a file named `clusterVector.jld2` being created, containing the structural data from the BHA run in the correct format to be used by the DACA for training. Simply copy this to the directory where the DACA is to be run, and rename the file to `explorationData.jld2`.