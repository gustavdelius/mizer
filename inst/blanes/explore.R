# Load parameters
params <- readRDS("inst/blanes/params.rds")

# Load catch distribution
catchdist <- readRDS("inst/blanes/catchdistribution.rds")

# Run tuning gadget
params <- tuneParams(params, catchdist = catchdist)
