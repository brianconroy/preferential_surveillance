# ------------------------------------- %
# PRISM variables and PCA
#
# Author: Ian Buller
# Date created: November, 8 2018
#
# ------------------------------------- %

# Setup -------------------------------------------------------------------
### Packages
library(prism) # prism data
library(sp) # spatial data manipulation
library(raster) # raster data manipulation
library(RStoolbox) # PCA of rasters
library(maps) # visualize geographical data
library(rgeos) # calculate buffers

### Set seed
set.seed(42) # the answer to life, the universe and everything

# Helper function for range transformation
range_scale <- function(reproj_rast){
  min_value <- min(na.omit(reproj_rast@data@values))
  max_value <- max(na.omit(reproj_rast@data@values))
  
  (reproj_rast - min_value) / (max_value - min_value)
}

crs_us <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Environmental Data ------------------------------------------------------

# Here, PRISM 30-Year Average Normals (can use bioclim or others)
# Download files
options(prism.path = "data/prism_raw")

prism_list     <- c("tmean", "tmax", "tmin", "ppt", "vpdmin", "vpdmax", "tdmean")
resolution_arg <- "4km"
annual_arg     <- TRUE
keepZip_arg    <- FALSE

# Load data
invisible(sapply(prism_list, FUN = get_prism_normals, 
                 resolution = resolution_arg, 
                 annual = annual_arg, 
                 keepZip = keepZip_arg)
          )

# Pre-processing Rasters ------------------------------------------------------
data_names <- ls_prism_data()$files
data_paths <- ls_prism_data(absPath = TRUE)$abs_path

# Convert to rasters
data_rasters         <- lapply(data_paths, raster)
names(data_rasters)  <- data_names

# Re-project Rasters
reproj_rasters       <- lapply(data_rasters, projectRaster, crs = crs(crs_us))

# Range transformation and stack
stand_reproj_rasters <- lapply(reproj_rasters, range_scale)
rasters_stand        <- stack(stand_reproj_rasters)

# PCA ---------------------------------------------------------------------

# Spatial PCA
pca1 <- RStoolbox::rasterPCA(rasters_stand)
summary(pca1$model) # PCA components
pca1$model$loadings # PCA loadings

# Extract Bands from PCA
pc1 <- pca1$map
pc1_b1 <- pc1[[1]] # PC1
pc1_b2 <- pc1[[2]] # PC2

# Mask scaled rasters by study area (window)
us <- raster::getData("GADM", country = "USA", level = 1)
ca <- us[match(toupper("California"),toupper(us$NAME_1)),]
# Extract outline in order to create buffer to capture all of PRISM
# region union kills the data frame so don't overwrite 'wus'
regs <- rgeos::gUnaryUnion(ca)
# takes way too long to plot without simplifying the polygons
regs <- rgeos::gSimplify(regs, 0.05, topologyPreserve = TRUE)
#plot(regs)
# Add 0.1 degree buffer to capture all of PRISM
ca_buffer <- rgeos::gBuffer(regs, width = 0.1, byid = TRUE) # same projection as crs_us

# Mask for California
mask_pc1 <- mask(pc1_b1, ca_buffer)
mask_pc2 <- mask(pc1_b2, ca_buffer)

ca_pc1 <- crop(pc1_b1, ca_buffer)
ca_pc1 <- mask(ca_pc1, ca_buffer)

ca_pc2 <- crop(pc1_b2, ca_buffer)
ca_pc2 <- mask(ca_pc2, ca_buffer)

plot(ca_pc1)
plot(ca_pc2)

ca_pcs <- stack(ca_pc1, ca_pc2)

writeRaster(ca_pcs,
            file = "data/prism_pcas_ca.grd", 
            bandorder = 'BIL',
            overwrite = TRUE)
