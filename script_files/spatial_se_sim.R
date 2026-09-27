#===============================================================================
# AEM 6850
# OLS, Conley and spatial error model standard errors under spatial dependence
#===============================================================================

#===============================================================================
# 1). Preliminary ---------
#===============================================================================

rm(list=ls())

# Load packages
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
packagelist <- c("spdep","plm","splm","RColorBrewer","maps","maptools","raster"
                 ,"data.table","geosphere","foreign","reshape","parallel","lfe","plyr","sf")
newpackages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(newpackages)) install.packages(newpackages)
sapply(packagelist, function(i) require(i, character.only=TRUE))
rm(newpackages, packagelist)

# Directories
dir <- list()
dir$root <- dirname(getwd())
dir$data <- paste(dir$root,"/data",sep="")
dir$output_figure <- paste(dir$root,"/output_figure",sep="")

# Misc
misc <- list()
misc$statefips <- c(23,33,50,44,25,9,11,36,34,10,24,51,54,39,42,21,17,18,55,26) 

# Functions
source("functions.R")

#===============================================================================
# 2). Downloading and cleaning the data -----
#===============================================================================

co <- map('county', plot=FALSE, fill=TRUE)
datareal <- data.frame(id=1:length(co$names), name= sapply(strsplit(co$names,":"), function(x) x[1])  )
co <- map2SpatialPolygons(co, datareal$id, proj4string = CRS(as.character(NA)))
co <- SpatialPolygonsDataFrame(co, datareal)
crs(co) <- CRS("+proj=longlat +datum=WGS84")
co@data$state <- sapply(strsplit(as.character(co@data[,2]),","), function(x) x[1])
co@data$county <- sapply(strsplit(as.character(co@data[,2]),","), function(x) x[2])
statevec <- sapply(strsplit(state.fips$polyname,":"),function(x) x[1])
co@data$statefips <- state.fips$fips[match(co@data$state,statevec)]
co@data$fips <- county.fips$fips[match(co@data[,2],county.fips$polyname)]
co <- co[co@data$statefips %in% misc$statefips,]

file <- list.files(dir$data, full.names=TRUE)
file <- file[grepl(".csv",file)]
dist <- read.csv(file)
dist <- dist[dist$HistoryFlag==1,]
dist$distfips <- dist$StateFips*10^3 + dist$DistrictCode
dist$fips <- dist$StateFips*10^3 + dist$CountyCode
co$distfips <- dist$distfips[match(co$fips, dist$fips)]

coords <- coordinates(co) 
co@data$lat <- coords[,2]
co@data$lon <- coords[,1]
n <- 958  # fixed number of counties
co <- co[1:n,]  # subset to first 958 counties

co <- st_as_sf(co)

# Spatial neighbors
sf_use_s2(FALSE)
nb <- knn2nb(knearneigh(coords[1:n,], k = 10), row.names = co$id)

# Spatial weights matrix
dlist <- nbdists(nb, coords[1:n,], longlat=T)  
idlist <- lapply(dlist, function(x) 1/x)
nbw <- nb2listw(nb, glist=idlist, style="W")
W_n <- listw2mat(nbw)

# Generate spatial correlation
years <- 2005:2009     # 5 years
t <- length(years)
rho <- 0.95
I_t <- diag(1,t)
sig2.u <- 25
I_n <- diag(1,n)	
M_n <- (I_n - rho * W_n)
Sigma <- sig2.u * (I_t %x% solve(t(M_n)%*%M_n))
Sigma1 <- Sigma[1:n,1:n]
vs <- svd(Sigma1) 
vcsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d))) 

# True regression coefficients
beta <- matrix(c(1,1), nrow=2, ncol=1)

# Create simulated dataset
data <- data.frame(
  id  = as.numeric(paste(co$id)), 
  fips= co$fips, 
  year= rep(years, each=n), 
  lat = co$lat, 
  lon = co$lon
)
data$distfips <- dist$distfips[match(data$fips, dist$fips)]

set.seed(123)
data$x2 <- mrnorm(t , rep(1,times=n)) # spatially dependent
data$x1 <- rnorm(t*n, sd=sd(data$x2)) # not spatially dependent
data$alpha <- rep(rnorm(n,mean=50,sd=sqrt(30)), times=t)

# Generate spatially-dependent errors
data$e <- mrnorm(t, rep(0,times=n))

#===============================================================================
# 3). Monte Carlo Simulation
#===============================================================================
nsims <- 750
outfile <- "out.RDS"

if(file.exists(outfile)) {
  out <- readRDS(outfile)
  message("Loaded existing simulation results from out.RDS")
} else {
  # Run MC simulation in parallel
  out <- mclapply(1:nsims, function(r) {
    
    # 1. Generate error and Y
    e <- mrnorm(t, rep(0, times=n)) 
    data$y <- as.matrix(data[,c("x1","x2")]) %*% beta + data$alpha + e   
    
    # 2. Run regression for OLS
    reg_ols  <-  plm(formula = y ~ x1 + x2, data = data, model="within", index = c("id","year"))
    
    # 2. Save coefficients and OLS SEs
    b_ols  <- coefficients(reg_ols)
    se_ols <- sqrt(diag(vcov(reg_ols)))
    names(se_ols) <- paste("se.", names(b_ols), sep="")
    
    # 3. Conley correction
    vcov_conley <- conley(reg_ols, idvec = data$id, timevec = data$year, 
                          latvec = data$lat, lonvec = data$lon, 
                          kernel = "bartlett", dist_cutoff = seq(0,500,100))
    se_conley <- sapply(vcov_conley[-1], function(x) sqrt(diag(x)))  # exclude OLS vcov
    # Save all SEs for each Conley cutoff 
    coef_names <- rownames(se_conley)
    cutoff_names <- paste0("cutoff", seq(0, 500, 100))
    colnames(se_conley) <- cutoff_names
    se_conley_vec <- as.vector(se_conley)
    names(se_conley_vec) <- as.vector(outer(coef_names, cutoff_names, paste, sep="_"))
    
    # 4. SEM
    reg_sem <- spml(formula = y ~ x1 + x2, data = data[,-2], model="within", index = c("id","year"),
                    effect = "individual", listw = nbw, lag = FALSE, spatial.error = "b")
    b_sem <- coefficients(reg_sem)
    se_sem <- sqrt(diag(vcov(reg_sem)))
    
    # Return results
    c(b_ols, se_ols, se_conley_vec, b_sem, se_sem)
    
  }, mc.cores = 6) # Use 6 cores
  
  out <- do.call("rbind", out)
  out <- as.data.frame(out)
  # Save the results to RDS
  saveRDS(out, file = outfile)
  message("Simulation results saved to out.RDS")
}

#===============================================================================
# 4). Plotting figure and extracting PNG
#===============================================================================

# Save as png with longer height than width
png(paste0(dir$output_figure, "/OLS_verses_SEM_Spatial_Dependence.png"), width=850, height=1400, res = 150)

# Rename output, columns 18-19 to SEM betas, and 21-22 to SEM SEs.
colnames(out)[18:19] <- c("b_sem.x1", "b_sem.x2")
colnames(out)[21:22] <- c("se_sem.x1", "se_sem.x2")

cutoffs <- seq(0, 500, 100) # Set up Conley cuttoffs for (0-500)
cols_red <- brewer.pal(length(cutoffs)+2, "Reds")[3:(length(cutoffs)+2)] # Set up color scale for Conley

par(mfrow=c(3,2), mar=c(5,4,2,1))


# ---- Top Panels Estimates of Betas 1 and 2 ----
ols_b1 <- density(out[,"x1"])
sem_b1 <- density(out[,"b_sem.x1"])
# Set x and y axis range for plotting beta1
y_max_beta1 <- max(ols_b1$y, sem_b1$y)
x_rng_beta1 <- range(c(ols_b1$x, sem_b1$x))
# Plotting beta 1
plot(ols_b1, col="red", lwd=2, main="", xlab="", ylim=c(0, y_max_beta1 * 1.1), xlim=x_rng_beta1, las=1)
lines(sem_b1, col="blue", lwd=2)
abline(v=beta[1], lty=2, lwd=1.5) # True beta 1
mtext("Estimates of beta1", side=1, line=3, cex=1.25)
legend("topleft", legend=c("OLS","SEM","True value"), 
       col=c("red","blue","black"), lty=c(1,1,2), cex=0.8, inset=c(0.02,0.02))

ols_b2 <- density(out[,"x2"])
sem_b2 <- density(out[,"b_sem.x2"])
# Set x and y axis range for plotting beta2
y_max_beta2 <- max(ols_b2$y, sem_b2$y)
x_rng_beta2 <- c(0.90, 1.10)
# Plotting beta 2
plot(ols_b2, col="red", lwd=2, main="", xlab="", ylim=c(0, y_max_beta2 * 1.1), xlim=x_rng_beta2, las=1, xaxt="n")
lines(sem_b2, col="blue", lwd=2)
abline(v=beta[2], lty=2, lwd=1.5) # True beta 2
# Custom x-axis labels
axis(side=1, at=seq(0.90, 1.10, by=0.05))
mtext("Estimates of beta2", side=1, line=3, cex=1.25)


# ---- Middle Panels SEM SEs ----
se1 <- density(out[,"se_sem.x1"])
true_sd_sem1 <- sd(out[,"b_sem.x1"]) # SD of SEM sampling distribution
# Force x-axis to start at 0 and go slightly past the max density value
x_max_se1 <- max(se1$x) * 1.5 
# Plotting SE beta 1 with custom x-axis
plot(se1, col="blue", lwd=2, main="", xlab="", xlim=c(0, 0.0065), ylim=c(0, max(se1$y)*1.1), las=1, xaxt="n")
abline(v=true_sd_sem1, lty=2, col="black") 
ticks <- seq(0, 0.006, by=0.001)
labels <- rep("", length(ticks))
labels[seq(1, length(ticks), by=2)] <- sprintf("%.3f", ticks[seq(1, length(ticks), by=2)])
axis(side=1, at=ticks, labels=labels)
mtext("Estimates of SE(beta1)", side=1, line=3, cex=1.25)

# Plotting SE beta 2
se2 <- density(out[,"se_sem.x2"])
true_sd_sem2 <- sd(out[,"b_sem.x2"])
x_max_se2 <- max(se2$x) * 1.5 
plot(se2, col="blue", lwd=2, main="", xlab="", xlim=c(0, x_max_se2), ylim=c(0, max(se2$y)*1.1), las=1)
abline(v=true_sd_sem2, lty=2, col="black")
mtext("Estimates of SE(beta2)", side=1, line=3, cex=1.25)


# ---- Bottom Panels Conley SEs ----
# Create function to handle Conley overlay plots and insert into plot_conley
plot_conley <- function(var_name, true_beta_col, xlim_vals=NULL) {
  # Identifying columns for variable's cutoffs
  col_names <- grep(paste0("^", var_name, ".*_cutoff"), names(out), value=TRUE)
  densities <- lapply(col_names, function(x) density(out[,x]))   # Calculate densities
  
  # Calculate Plot Limits
  ymax <- max(sapply(densities, function(d) max(d$y)))   # y-max, is the highest peak among all cutoff lines
  xmax <- max(sapply(densities, function(d) max(d$x)))   # x-max, is the furthest point to the right
  
  # Plot, first cutoff
  if(is.null(xlim_vals)){
    plot(densities[[1]], type="l", col=cols_red[1], lwd=1.5,
         ylim=c(0, ymax*1.1), xlim=c(0, xmax*1.1),
         main="", xlab="", las=1)
  } else {
    plot(densities[[1]], type="l", col=cols_red[1], lwd=1.5,
         ylim=c(0, ymax*1.1), xlim=xlim_vals,
         main="", xlab="", las=1)
  }
  
  for(i in 2:length(densities)){
    lines(densities[[i]], col=cols_red[i], lwd=1.5)
  }   # Add remaining lines
  abline(v=sd(out[,true_beta_col]), lty=2, col="black")   # Add true SD
  # Adjust x-axis lables
  if(var_name=="x1"){
    mtext("Estimates of SE(beta1)", side=1, line=3, cex=1.25)
  } else if(var_name=="x2"){
    mtext("Estimates of SE(beta2)", side=1, line=3, cex=1.25)
  }
}
# Plotting Conley beta 1 and 2 with specified axes
plot_conley("x1", "x1", xlim_vals=c(0, 0.018)) # bottom left
legend("topleft", title="Conley cutoff (mi):", 
       legend=cutoffs, col=cols_red, lty=1, cex=0.7, ncol=2, inset=c(0.02,0.02))
plot_conley("x2", "x2", xlim_vals=c(0, 0.10)) # bottom right
dev.off()

