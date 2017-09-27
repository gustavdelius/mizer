params_data2 <- read.csv("./vignettes/NS_species_params.csv")
params_data2

load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]

tf <- t(Fmat)
tf
rownames(tf)[19:29]

dim(tf)

avfromFmat <- colMeans(tf[19:29,])

baseF <-c(1.0399246, 0.8685833, 1.4873488, 0.4000000, 0.7665488, 0.4000000,
          0.8171929, 1.1308975, 0.6881135, 1.1435399, 0.9045996, 0.9882556)

# agrees for Northen pout and Herriing

#############

f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
f_history_means <- as.numeric(colMeans((f_history)[19:29,]))


######

mikes_species_list <- c("Sprat",
                        "Sandeel",
                        "N.pout",
                        "Dab",
                        "Herring",
                        "Gurnard",
                        "Sole",
                        "Whiting",
                        "Plaice",
                        "Haddock",
                        "Saithe",
                        "Cod"
)

getmsindex <- function(s){
  matcher <- match(mikes_species_list,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}
getmsindex("Cod")

getnsindex <- function(s){
  matcher <- match(params_data$species,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}

baseF[1]

baseFpermtoNS <- 1:12
for (i in (1:12)){
  baseFpermtoNS[i] <- baseF[getnsindex(mikes_species_list[i])]
}

avfromFmattoNS <- 1:12
for (i in (1:12)){
  avfromFmattoNS[i] <- avfromFmat[getnsindex(mikes_species_list[i])]
}


round(f_history_means,2)
round(baseFpermtoNS,2)
round(avfromFmattoNS,2)

dim(f_history)
f_history[,1]
Fmat[1,]
