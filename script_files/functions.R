# Functions
# Another functions that makes it easy to do multivariate normal draws for a given covariance matrix
mrnorm <- function(num.vals, mu.vec) {
  k        <- ncol(Sigma1)
  ans.mat  <- sweep(matrix(rnorm(num.vals*k), nrow = num.vals) %*% vcsqrt,2,mu.vec,"+") 
  c(t(ans.mat))
} 
# Spatial HAC (Conley)
conley <- function(reg, idvec, timevec, latvec, lonvec, kernel = "bartlett",dist_cutoff = 500, ncores=1) {
  # Syntax: conley(reg, data$id, data$year, data$lat, data$lon, kernel = "bartlett",dist_cutoff = seq(0,500,100), ncores=1)
  # reg: a plm or felm object
  # idvec: a vector of location IDs (n x t)
  # timevec: a vector of time IDs (n x t)
  # latvec: a vector of latitute (n x t) 
  # lonvec: a vector of longitude (n x t) 
  # kernel: one of "unfiform", "bartlett" and "bartlettprod"
  # dist_cutoff: a vector of distance cutoffs (miles)
  
  # Load packages
  wants <- c("data.table","geosphere","foreign","reshape","parallel","plm","lfe","plyr")
  has   <- wants %in% rownames(installed.packages())
  if(any(!has)) install.packages(wants[!has])
  sapply(wants, function(i) require(i, character.only=TRUE))
  
  # Internal functions (note:distance is in meters)
  weightfun <- list(uniform=function(dist, cut) {dist<=cut}, bartlett=function(dist, cut) {w <- dist<=cut ; (1-dist/(cut+1)) * w } )
  Fac2Num <- function(x) {as.numeric(as.character(x))}
  iterateObs   <-function(y1,e1,X1,fordist,coefficients,cutoff=500) {
    # Distances and kernel weight
    distances2 <- as.matrix(dist(fordist)) * 69 # in miles
    weights2  <- apply(distances2, 1, function(x) weightfun[[kernel]](dist=x, cut=cutoff)) 
    E1 <- t(t(e1)) %*% t(e1) 
    XeeXhs <- t(X1) %*% (E1 * weights2) %*% X1
  } 
  
  if(class(reg)[1] == "felm") {
    Xvars <- rownames(reg$coefficients)
    dat1 <- data.frame(reg$cY, reg$cX, e=c(reg$residuals))
    dat2 <- data.frame(id=idvec, time=timevec, lat=latvec, lon=lonvec) 
    dat2 <- dat2[order(dat2$time,dat2$id),]
    dat <- cbind(dat1,dat2)  
    olsvcov <-  reg$vcv
  } else if (class(reg)[1] == "plm") {
    Xvars <- names(reg$coefficients)
    dat2 <- data.frame(id=idvec, time=timevec, lat=latvec, lon=lonvec) 
    dat2 <- dat2[order(dat2$id,dat2$time),]
    temp <- data.frame(id=idvec, time=timevec)
    temp <- temp[order(temp$id, temp$time),]
    temp <- as.data.frame(cbind(temp, reg$model))
    temp <- ddply(temp, .(id), numcolwise(scale, scale = FALSE))
    temp <- temp[order(temp$id, temp$time),]
    dat1 <- data.frame(temp[,-(1:2)], e=c(reg$residuals))
    dat <- cbind(dat1,dat2)
    dat <- dat[order(dat$time, dat$id),]
    names(dat)[2:(length(Xvars)+1)] <-  Xvars
    olsvcov <-  reg$vcov
  } else {
    message("Incompatible panel model class. Use felm or plm objects.")
    break
  }
  
  n <- nrow(dat)
  k <- length(Xvars)
  
  # Correct for spatial correlation:
  timeUnique <- unique(dat$time)
  Ntime <- length(timeUnique)
  coef <- coefficients(reg)
  
  # Loop over cutoffs
  out <- lapply(dist_cutoff, function(d) {
    
    print(paste(round(d,1),"miles"))
    
    # Compute
    if (ncores>1) {
      XeeXhs <- mclapply(timeUnique, mc.cores=ncores, FUN=function(t) iterateObs(dat$y[dat$time==t],dat$e[dat$time==t],as.matrix(dat[dat$time==t,Xvars,drop=F]),as.matrix(dat[dat$time==t,c("lon","lat")]),coef,cutoff=d))
    } else {
      XeeXhs <- lapply(timeUnique, function(t) iterateObs(dat$y[dat$time==t],dat$e[dat$time==t],as.matrix(dat[dat$time==t,Xvars,drop=F]),as.matrix(dat[dat$time==t,c("lon","lat")]),coef,cutoff=d))
    }
    XeeX <- Reduce("+",  XeeXhs)
    
    # Generate VCE for only cross-sectional spatial correlation:
    X <- as.matrix(dat[, Xvars])
    invXX <- solve(t(X) %*% X) * n
    V <- invXX %*% (XeeX / n) %*% invXX / n
  })
  
  # Export
  return_list <- c( list(olsvcov), out)
  names(return_list) <- c(-1, dist_cutoff)
  return(return_list)
}
