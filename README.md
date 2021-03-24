
# Preferential Surveillance

This repository contains code to estimate disease risk surfaces
from preferentially sampled disease surveillance data. It includes simulation studies
comparing performance of the proposed method against various benchmarks
(Bayesian Additive Regression Trees, Spatial Poisson Regression, Poisson Regression)
in estimating disease risk from preferentially sampled data. Also included is an analysis 
of real surveillance data monitoring plague (Yersinia pestis) within the sciurid population of California. 

# Organization

This repository can be installed like a package for ease of access to
function documentation, processed datasets, and dependency management.
However it also contains a number of "project like" elements including
simulation and analysis scripts. 

Simulation scripts are organized as follows:

* **simulation_data.R**: generates simulated disease surveillance datasets under differing
levels of preferential sampling.
* **simulation.R**: Fits the proposed preferential sampling method along with benchmarks to the 
datasets simulated by simulation_data.R. 
* **simulation_analysis.R**: Compares performance of the models fit by **simulation.R**.

Plague analysis scripts are structured as:

* **data_analysis_fits.R**: Fits the proposed preferential sampling method along with benchmarks to the 
processed plague surveillance dataset. Additionally, performs spatial downscaling to obtain high resolution
risk maps.
* **data_analysis_summary.R**: Summarizes, visualizes, and compares risk maps obtained from the proposed
method and benchmarks over theplague surveillance dataset.

Remaining files are organized as in any other package:

* **R/**: holds all helper functions
* **data-raw**: contains scripts to generate the processed datasets of the package
* **man/**: documentation

# Installation

Dependencies of the repository and proccessed datasets can be installed and made accessible by
first cloning the repository, navigating to the root of the repository in the R console, and 
installing the repository locally via:

```{r}
devtools::install()
```

Successful installation can be checked by loading the package and inspecting its datasets.

```{r}
library(preferentialSurveillance)

# first and second principal components of the PRISM climatic dataset
plot(prism_pca)
```