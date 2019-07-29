params <- readRDS("inst/humboldt/params2607.rds")

dx <- log10(params@w[2]/params@w[1])

ppmr <- read_csv("inst/humboldt/Data_PPMR.csv") %>% 
    mutate(Species = str_to_title(Species))

species <- "Sardine"
stomach <- filter(ppmr, Species == species)
species_params <- params@species_params[params@species_params$species == species, ]
pl <- fitKernel(stomach,
                species_params = species_params,
                dx = dx, lambda = params@lambda)

ppmr <-
    read_tsv("inst/ppmr/Predator_and_prey_body_sizes_in_marine_food.csv",
             na = c("", "n/a"),
             guess_max = 10000) %>% 
    select(Species = `Predator common name`,
           wpredator = `SI predator mass`,
           wprey = `SI prey mass`)

pl <- fitKernel(filter(ppmr, Species == "Swordfish"),
                dx = dx, lambda = 2)

