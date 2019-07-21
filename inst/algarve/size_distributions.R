#' Load the data
library(tidyverse)
Individual_fish_sizes <- read_csv("Individual fish sizes.csv")
params <- readRDS("params_102_calibrated.rds")

#' We massage the data frame with the sizes a bit
sizes <- Individual_fish_sizes %>%
    rename(Species = `common name`) %>% 
    mutate(Species = recode(Species, "Dourada" = "Dourada "))%>% 
    group_by(Species) 

#' Note that for some of the species some observed sizes are larger than 
#' the value of w_inf, which should probably be fixed.
max_w <- sizes %>% 
    summarise(max_w = max(weight)) %>% 
    mutate(w_inf = params@species_params[Species, "w_inf"])
max_w

#' Now let us look only at the species with more than 10 observations
sizes <- sizes %>% 
    filter(n() > 10)
unique(sizes$Species)

#' Next we want to plot the observed size distribution together with the
#' model size distribution. To do that with ggplot we first need to create
#' a data frame with the model distributions.
modeldens <- plyr::ddply(sizes, "Species", function(df) {
    sp <- unique(df$Species)
    idx_min <- sum(params@w < min(df$weight))
    idx_max <- sum(params@initial_n[sp, ] > 0) + 5
    select <- idx_min:idx_max
    x <- log(params@w[select])
    dx <- x[2] - x[1]
    density <- params@initial_n[sp, select] * params@w[select]
    density <- density / sum(density) / dx
    data.frame(x, density)
})

#' Now we can plot
ggplot(sizes, aes(log(weight))) +
    geom_density() +
    facet_wrap(~`Species`, scales = "free") +
    geom_line(aes(x, density), data = modeldens,
              colour = "blue")


