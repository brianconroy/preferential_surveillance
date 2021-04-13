
# Preferential Surveillance

This repository contains code to estimate disease risk surfaces
from preferentially sampled disease surveillance data. It includes simulation studies
comparing performance of the proposed method against various benchmarks
(Bayesian Additive Regression Trees, spatial Poisson regression, Poisson regression). 
Also included is an analysis of real surveillance data monitoring plague (*Yersinia pestis*) 
within the rodent population of California. 

This repository can be installed like a package for ease of access to
function documentation, processed datasets, and dependency management.
However it also contains a number of "project like" elements including
simulation and analysis scripts. 

Simulation scripts are organized as follows:

* **simulation_data.R**: generates simulated disease surveillance datasets under differing
levels of preferential sampling.
* **simulation.R**: Fits the proposed preferential sampling method along with benchmarks to the 
datasets simulated by simulation_data.R. 
* **simulation_analysis.R**: Compares performance of the models fit by simulation.R.

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

**Assumptions**:

* R version 3.4.3 - 4.0.5 has been installed.
* The devtools package has been installed.

```{r}
install.packages("devtools")
```

Dependencies of the repository and proccessed datasets can be installed and made accessible by
first cloning the repository, navigating to the root of the repository in the R console, and 
installing the repository locally via:

```{r}
# install package dependencies
devtools::install_deps()

# install package
devtools::install()
```

Successful installation can be checked by loading the package and inspecting its datasets.

```{r}
library(preferentialSurveillance)

# first and second principal components of the PRISM climatic dataset
plot(prism_pca)
```

# Data Access

**Caveat**. The scripts for data analysis in this repo consume a plague surveillance
dataset of aggregate (1983-2015) counts of plague positive and
negative rodents. This dataset was collected and is maintained by the California Department of Public Health - 
Vector Borne Disease Section.
Due to department policy the data are **not** included in the repository, but are available upon request.
