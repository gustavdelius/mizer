library(mizer)
library(plyr)




source("./R/project_methods.R")
source("./R/project.R")
paramsConst <- set_trait_model(no_sp = 8, min_w_inf = 10, max_w_inf = 1e5)
simConst <- project(paramsConst, t_max=50, effort = 0)
plot(simConst)
####################################################################################

epsi <- 0.001
source("./R/Experimental Code/project_methodsmodPREYSWITCH.R")
source("./R/Experimental Code/projectmodPREYSWITCH.R")
paramsConst <- set_trait_model(no_sp = 8, min_w_inf = 10, max_w_inf = 1e5)
simConst <- project(paramsConst, t_max=50, effort = 0)
plot(simConst)

###################################################################### free versions
source("./R/project_methods.R")
source("./R/project.R")
paramsConst <- set_trait_model(no_sp = 8, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simConst <- project(paramsConst, t_max=50, effort = 0)
plot(simConst)
####################################################################################

epsi <- 0.1
source("./R/Experimental Code/project_methodsmodPREYSWITCH.R")
source("./R/Experimental Code/projectmodPREYSWITCH.R")
paramsConst <- set_trait_model(no_sp = 8, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simConst <- project(paramsConst, t_max=50, effort = 0)
plot(simConst)
#############



####
paramsFree <- set_trait_model(no_sp = 2, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simFree <- project(paramsFree, t_max=75, effort = 0)
plot(simFree)

####################

paramsConst2 <- set_trait_model(no_sp = 3, min_w_inf = 10, max_w_inf = 1e5)
simConst2 <- project(paramsConst2, t_max=750, effort = 0)
plot(simConst2)


paramsFree2 <- set_trait_model(no_sp = 2, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simFree2 <- project(paramsFree2, t_max=75, effort = 0)
plot(simFree2)

###

paramsConst3 <- set_trait_model(no_sp = 3, min_w_inf = 10, max_w_inf = 1e5)
simConst3 <- project(paramsConst3, t_max=750, effort = 0)
plot(simConst3)

paramsFree3 <- set_trait_model(no_sp = 3, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simFree3 <- project(paramsFree3, t_max=75, effort = 0)
plot(simFree3)


####
