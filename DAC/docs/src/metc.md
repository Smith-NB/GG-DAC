# Metropolis Criteria

```@meta
CurrentModule = DAC
```

Metropolic Criteria types must implement the method `getAcceptanceBoolean(MetC::EnergyMetC, oldCluster::Cluster, newCluster::Cluster)` which returns a `Tuple` with the acceptance or rejection of the given step (i.e. `true` or `false`), and a `String` for debugging that denotes the chance to accept the hop (if any) and any other relavent information.

# Energy Metropolis Criterion

```@docs
EnergyMetC
```

# HISTO Metropolis Criterion

```@docs
HISTOMetC
HISTOMetC(kT::Float64, w::Float64, delta::Float64, resetPeriod::Float64, clusterVector::ClusterVector, waitTime::Int64, io::Tuple{IO, Channel})
```

# GMM Metropolis Criterion
```@docs
GMMMetC
GMMMetC(gaussian::GMM, gaussianCluster::Int64, pca::PCA, mode::Symbol, useExplorationDataOnly::Bool, kT::Float64, io::Tuple{IO, Channel})
GMMMetC(gaussian::GMM, gaussianCluster::Int64, pca::PCA, mode::Symbol, useExplorationDataOnly::Bool, kT::Float64, classes::normalCNAProfile, io::Tuple{IO, Channel})```
