
#######################################################################################

# This script processes the raw plague surveillance dataset (1983-2015) from 
# the California Department of Public Health - Vector Borne Disease Section
# The output dataset contains aggregated counts of cases (plague positive 
# sciurids), controls (plague negative sciurids). Specifically, the output
# is a list of the following lists:

# NOTE: The raw surveillance dataset is not publicly accessible due to 
# privacy concerns (relating to property values and the exact locations 
# of plague positive rodents...). This script will only run if the raw
# dataset is locally accessible.

# locs: a list describing the locations at which rodents were sampled for plague
#       over the (rasterized) state of California. Contains the following values:
#         * raster: a raster of the state of CA at the desire resolution of analysis (426 km^2),
#                   whose values are binary indicators for whether the cell was observed by the 
#                   surveillance system
#         * cells: cell ids of the raster which have been sampled by the surveillance system
#         * status: vector of binary values (0, 1) indicating whether each cell in the raster
#                   has been observed
#         * coords: matrix of lat, lon coordinates of the center points of raster cells that 
#                   have been observed
#         * ids: ids of the raster which have been sampled by the surveillance system where 
#                ids are identified after excluding NA values of the raster

# case.data: a list describing case counts (plague positive) over the study area. Contains 
#           the following values:
#         * y: vector of counts of plague positive specimen recovered over the observed 
#              raster cells
#         * x.standardised: design matrix containing an intercept and two columns, the 
#                           first and second principal components of the PRISM dataset
#                           after range and scale standardisation
#         * x: design matrix containing an intercept and two columns, the 
#              first and second principal components of the PRISM dataset
#              (not standardised)
#         * p: equal to the number of covariates plus intercept

# ctrl.data: a list describing control counts (plague negative) over the study area. 
#            Contains the following values:
#         * y: vector of counts of plague negative specimen recovered over the observed 
#              raster cells.
#         * x.standardised: design matrix containing an intercept and two columns, the 
#                           first and second principal components of the PRISM dataset
#                           after range and scale standardisation
#         * x: design matrix containing an intercept and two columns, the 
#              first and second principal components of the PRISM dataset
#              (not standardised)
#         * p: equal to the number of covariates plus intercept

#######################################################################################

library(raster)
library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
library(jsonlite)
library(fields)
#library(preferentialSurveillance)
sourceDirectory('R/')

## Define Global Variables ------------------------------------------------------------

# Load the raw surveillance dataset. This dataset is not publicly
# accessible due to privacy concerns (relating to property values
# and the exact locations of plague positive rodents...)
SRC <- "/Users/brianconroy/Documents/research/cdph/data/CDPH_scurid_updated_full.csv"

## Load Data --------------------------------------------------------------------------

# Load and aggregate the PRISM principal components
caPr <- prism_pca
caPr.disc <- aggregate(caPr, fact=5)

# Check the number of grid cells
N <- n_values(caPr.disc[[1]])
print(N)

# Load the raw surveillance dataset
rodents <- read.csv(SRC, header=T, sep=",")

# Create a raster outlining the aggregated studdy area
loc.disc <- caPr.disc[[1]]

## Create Processed Dataset -----------------------------------------------------------

# Create the dataset of locations, case and control counts
plague_data <- assemble_data(rodents, caPr.disc)

# Check the raster indicating location of sample sites
plot(plague_data$loc$raster)

## Save Output ------------------------------------------------------------------------

# Save the data to be loading when the package loads
usethis::use_data(plague_data)
