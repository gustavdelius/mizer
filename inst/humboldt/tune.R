catch <- readRDS("~/Git/mizer/inst/humboldt/catch.rds")
humboldt_params <- readRDS("~/Git/mizer/inst/humboldt/humboldt_params.rds")

tuneParams(humboldt_params, catch = catch)


