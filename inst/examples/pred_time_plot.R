
######fn from http://www.phaget4.org/R/image_matrix.html
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
    min <- min(x)
    max <- max(x)
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    title <-c()
    # check for additional function arguments
    if( length(list(...)) ){
        Lst <- list(...)
        if( !is.null(Lst$zlim) ){
            min <- Lst$zlim[1]
            max <- Lst$zlim[2]
        }
        if( !is.null(Lst$yLabels) ){
            yLabels <- c(Lst$yLabels)
        }
        if( !is.null(Lst$xLabels) ){
            xLabels <- c(Lst$xLabels)
        }
        if( !is.null(Lst$title) ){
            title <- Lst$title
        }
    }
    # check for null values
    if( is.null(xLabels) ){
        xLabels <- c(1:ncol(x))
    }
    if( is.null(yLabels) ){
        yLabels <- c(1:nrow(x))
    }
    
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
    
    # Red and green range from 0 to 1 while Blue ranges from 1 to 0
    ColorRamp <- rgb( seq(0,1,length=256),  # Red
                      seq(0,1,length=256),  # Green
                      seq(1,0,length=256))  # Blue
    ColorLevels <- seq(min, max, length=length(ColorRamp))
    
    # Reverse Y axis
    reverse <- nrow(x) : 1
    yLabels <- yLabels[reverse]
    x <- x[reverse,]
    
    # Data Map
    par(mar = c(3,5,2.5,2))
    image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
          ylab="", axes=FALSE, zlim=c(min,max))
    if( !is.null(title) ){
        title(main=title)
    }
    axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
    axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
         cex.axis=0.7)
    
    # Color Scale
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    
    layout(1)
}
# ----- END plot function ----- #

MM <- matrix(c(1,2,3,4),nrow = 2,ncol=2)
myImagePlot(MM)

#to do

# get dynamics from sim

# extract dominant predators

# form them in a matrix of weight vs time, according to numerical index

# plot the result using 
#myImagePlot(MM)

#put sensible colors into the humboldt system

# do a custom color version of thid code, and a version for `most beneficial prey



######### another copy of the preamble from humboldt.R

params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))
effort <- 1.4

no_sp <- dim(params_data)[1]

l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)

p <- setBackground(
    set_scaling_model(min_w_pp = 1e-12,
                      no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
                      min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
                      lambda = 2.12,
                      knife_edge_size = Inf)
)


all_efforts <- c(0, 1.4, 1.1)
names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
effort <- all_efforts[1:2]
for (i in (1:no_sp)) {
    if (params_data$species[i] == "Anchovy") {
        effort <- c(effort, all_efforts[3])
    }
    a_m <- params_data$a2[i]
    b_m <- params_data$b2[i]
    L_inf_m <- params_data$Linf[i]
    L_mat <- params_data$Lmat[i]
    if (params_data$species[i] == "Anchovy") {
        gear <- "sigmoid_gear_Anchovy"
    } else {
        gear <- "sigmoid_gear"
    }
    species_params <- data.frame(
        species = as.character(params_data$species[i]),
        w_min = params_data$Wegg[i],
        w_inf = params_data$w_inf[i],
        w_mat = params_data$w_mat[i],
        beta = params_data$beta[i],
        sigma = log(params_data$sigma[i]),
        z0 = 0,
        alpha = 0.6,
        erepro = 0.1, # unknown, determined later
        sel_func = "sigmoid_length",
        gear = gear,
        l25 = l25[i],
        l50 = l50[i],
        k = 0,
        k_vb = params_data$k_vb[i],
        a = a_m,
        b = b_m
    )
    
    p <- addSpecies(p, species_params, effort = effort, rfac=Inf)
}

# Run to steady state
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-3)

###################



# 
# mu_pairwise <- function(sim,i,j){
#     
#     # for a given time t, let n be the current state, and let m be like n, 
#     # but with all species removed except i, now we want 
#     
#     
#     for (t in (1:dim(sim@n)[1])){
#         
#     }
#     
#     getPredRate(params,m,b_pp)[j,]
#         
#     
# }
# 

############################ make random initial condition


nn <- p@initial_n

for (i in 1:no_sp){
    nn[i, ] <- runif(length(nn[i, ]))
}
t_max <- 15
t_save <- 0.1
sim <- project(p, t_max = t_max, t_save = t_save, effort = effort, initial_n = nn)
plot(sim)


no_sp <- dim(sim@n)[2]



#########################

j <- 2
                        
      # time by weight                  
MM <- matrix(0,nrow = dim(sim@n)[1], ncol = dim(sim@n)[3])

for (t in (1:dim(sim@n)[1])){
    n <- sim@n[t,,]
    mumu <- sim@params@psi
    for (i in 1:no_sp){
        m <- n
        m[,] <- 0
        m[i,] <- n[i,]
        full_pr <- getPredRate(sim@params,m,sim@n_pp[t,])
        mumu[i,] <- full_pr[j,(dim(full_pr)[2]-length(sim@params@w)+1):dim(full_pr)[2]]
    }
    for (w in 1:dim(sim@n)[3]){
        MM[t,w] <- which.max(mumu[,w])
    }
    
}

myImagePlot(MM)


MM[]
                        
####################
#For a given species, I guess we can show how the dominant predators change over time by visualizing the colour changing on the weight line using a 2D diagram similar to the space time diagrams used for 1D cellular automata
#The diagram would be a grid where cell (x,y) is given the colour corresponding to the dominant predator species on weight y at time x
#

            `