params <- readRDS("inst/blanes/params.rds")

library(readr)
size_frequency <- read_csv("inst/blanes/size_frequency.csv")
size_frequency$cm <- 3:56

library(reshape2)
catchdistribution <- melt(size_frequency, id.vars = "cm", na.rm = TRUE)
levels(catchdistribution$variable) <- 
    c("Angler fish", "Hake", "Blue whiting", "Red mullet", "Striped red mullet",
      "Horse mackerel", "Poor cod", "Shortfin squid", "Gurnards", "Starfish")
names(catchdistribution) <- c("length", "species", "catch")
View(catchdistribution)

saveRDS(catchdistribution, "inst/blanes/catchdistribution.rds")

tuneParams(params, catchdist = catchdistribution)
