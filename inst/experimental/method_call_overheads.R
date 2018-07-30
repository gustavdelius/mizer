# The purpose of this speed comparison is to determine the overhead arising
# from the method calls during the project loop.

# The result shows that the overhead is negligible

library(microbenchmark)
library(mizer)
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)

microbenchmark(project(params, effort=1, t_max=200, fast = FALSE),
               project(params, effort=1, t_max=200, fast = TRUE),
               times = 10)
