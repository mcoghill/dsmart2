# This will demonstrate the capability of the terra implementation of dsmart, 
# as opposed to the original method of using the raster package.

# This example will create a run from the example data using the caret package,
# specifically using "ranger". First, a bulk load of all of the packages for
# this run is necessary:

ls <- c("tidyverse", "terra", "stars", "raster", "caret", "ranger", "MLmetrics")
new.packages <- ls[!(ls %in% installed.packages()[, "Package"])]
if(length(new.packages))
  install.packages(new.packages)

# Make sure terra package is up to date! This may take a moment
if(compareVersion(as.character(packageVersion("terra")), "1.0-7") < 0)
  install.packages("terra")
invisible(lapply(ls, library, character.only = TRUE)[0])
rm(ls, new.packages)

# Load custom DSMART files
source("./R/dsmart.R")

# Load DSMART example data
data(dalrymple_composition, package = "rdsmart")
data(dalrymple_observations, package = "rdsmart")
data(dalrymple_covariates, package = "rdsmart")
data(dalrymple_polygons, package = "rdsmart")

# Convert example data to work within new DSMART implementation
dalrymple_covariates <- terra::rast(dalrymple_covariates)
dalrymple_polygons <- terra::vect(dalrymple_polygons)

# Run dsmart with extra observations
outputdir <- file.path(getwd(), "dsmart2")
dir.create(outputdir, showWarnings = FALSE)
dsmart2 <- dsmart(
  covariates = dalrymple_covariates, 
  polygons = dalrymple_polygons, 
  composition = dalrymple_composition,
  observations = dalrymple_observations,
  rate = 15, reals = 10, outputdir = outputdir)
# That took slightly over 2 minutes
# Now, let's see how long it takes using the default package
# Remove functions so package doesn't get confused
rm(disaggregate, dsmart, summarise, .allocate, .getStratifiedVirtualSamples, .getVirtualSamples, .observations)

# Load package and associated data
library(rdsmart)
data(dalrymple_composition, package = "rdsmart")
data(dalrymple_observations, package = "rdsmart")
data(dalrymple_covariates, package = "rdsmart")
data(dalrymple_polygons, package = "rdsmart")

# Create separate output directory and run
outputdir <- file.path(getwd(), "dsmart_default")
dir.create(outputdir, showWarnings = FALSE)

dsmart_default <- rdsmart::dsmart(
  covariates = dalrymple_covariates, 
  polygons = dalrymple_polygons, 
  composition = dalrymple_composition,
  observations = dalrymple_observations,
  rate = 15, reals = 10, 
  outputdir = outputdir, 
  cpus = parallel::detectCores()
)
# That took 6 minutes, 3 times longer than dsmart2


# debugging
covariates = dalrymple_covariates
polygons = dalrymple_polygons
composition = dalrymple_composition
observations = dalrymple_observations
rate = 15 
reals = 10
method.sample = "by_polygon"
method.allocate = "weighted"
method.model = NULL
args.model = NULL
strata = NULL
nprob = 3
outputdir = file.path(getwd(), "dsmart2")
stub = NULL
factors = NULL
type = "raw"
tilesize = 500
outDir = file.path(outputdir, "tiles")
mask = NULL
