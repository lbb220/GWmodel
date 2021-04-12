# NB - need to reference the following journal papers for all this GWPCA work in GWmodel pdf manual:

#   	Harris P, Brunsdon C, Charlton M (2011)
#     Geographically weighted principal components analysis.
#     International Journal of Geographical Information Science 25:1717-1736

#     Harris P, Juggins S, Clarke A, Brunsdon C, Charlton M (2013)
#     Investigating spatial structure in an ecological data set using geographically weighted principal components analysis.
#     In review

#     Harris P, Brunsdon C, Charlton M, Juggins S, Clarke A (2013)
#     Multivariate spatial outlier detection using robust geographically weighted techniques.
#     In review

#     Harris P, Brunsdon C, Charlton M Juggins S, Clarke A (2013)
#     Robust geographically weighted principal components analsysis and its use in sample re-design.
#     In review
#################################################################################################################  
# 1. Basic or robust GWPCA functions ############################################################################
#################################################################################################################
##Author: PH
##Edited by BL
# Basic WPCA with bi-square weighting... 
wpca <- function(x,wt,...) {
	local.center <- function(x,wt)  
		sweep(x,2,colSums(sweep(x,1,wt,'*'))/sum(wt))
	svd(sweep(local.center(x,wt),1,sqrt(wt),'*'),...)}
  	
# Robust SVD function...
# This is done so that a direct replacement can be found for our original attempt using robustSvd() in pcaMethods R package (bioconductor)
# NB the true variances (d) can be returned unlike in the basic GWPCA case    i.e. d=(pc$sdev)^2
# However the squared term is removed in order to match the basic GWPCA case - this may need tidying up...
robustSvd <- function(x) {
  pc <- princomp(covmat=covMcd(x,alpha=3/4)$cov) # observe aplha set at 0.75
  return(list(v=pc$loadings,d=pc$sdev))}
# Robust WPCA with bi-square weighting (note using medians as a robust estimate of location)....
# First doing the medians...
wt.median <- function(x,wt) {
 wt.median.1 <- function(x,wt)	{
	ox <- order(x)
	wox <- cumsum(wt[ox])
	posn <- which.min(abs(wox - 0.5))
	return(x[ox][posn]) }
 return(apply(x,2,wt.median.1,wt))}


# Then the weighted PCA... 
rwpca <- function(x,wt,nu=0,nv=2) { 
  mids <- sweep(x,2,wt.median(x,wt))
  res <- robustSvd(sweep(mids,1,wt,'*'))
  res$v <- res$v[,1:nv]
  return(res)}



# Main GWPCA function
# x is the data matrix,
# loc is a 2-column coordinate matrix
# bw is the bandwidth in fixed or adaptive form,
# k is the number of retained components
# eloc are locations to evaluate the loadings,  if different from loc - eloc is also a 2 column matrix
# pcafun is the weighted PCA function which can be in a basic or robust form (from before)
# This returns a list with <loadings> for each eloc; d for each eloc; and the bandwidth used	
#gwpca <- function(x,loc,bw,k=2,eloc=loc,pcafun=wpca,...) 
gwpca <- function (data, elocat, vars, k = 2, robust = FALSE, kernel = "bisquare",
                  adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, cv = T, scores=F,
                  dMat)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    polygons <- NULL
    if(is(data, "SpatialPolygonsDataFrame"))
       polygons <- polygons(data)
  }
  else
     stop("Given data must be a Spatial*DataFrame or data.frame object")
  
  if (missing(elocat))
  {
  	ep.given <- FALSE
    elocat<-coordinates(data)
  }
  else 
  {
    ep.given <- T
    if (is(elocat, "Spatial"))
    {
       espdf<-elocat
       elocat<-coordinates(espdf)
    }
    else if (is.numeric(elocat) && dim(elocat)[2] == 2)
       elocat<-elocat
    else
      {
        warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
        elocat<-dp.locat
      }
  }
  data <- as(data, "data.frame")
  dp.n<-nrow(data)
  ep.n<-nrow(elocat)
  if (missing(dMat))
  {
    DM.given<-F
    DM1.given<-F
    if(dp.n + ep.n <= 10000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=elocat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    DM1.given<-T 
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=ep.n)
       stop("Dimensions of dMat are not correct")
  }
  if(missing(vars))
    stop("Variables input error")
  if(missing(bw)||bw<=0)
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  ##########Extract data
  col.nm<-colnames(data)
  var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if(length(var.idx)==0) stop("Variables input doesn't match with data")
  x<-data[,var.idx]
  x<-as.matrix(x)
  var.nms<-colnames(x)
  var.n<-ncol(x)
  if(len.var > var.n)
     warning("Invalid variables have been specified, please check them again!")
  #####PCA
  pca.res <- princomp(x, cor=T, scores=scores)
  ############################# WPCA
	#The local loading for the principle components
  w <- array(0,c(ep.n,var.n,k))
  #Return the local scores
  if(scores)
    gwpca.scores <- list()
  else
    gwpca.scores <- NULL
  #
	d <- matrix(0,ep.n,var.n)
	# Add this bit please ############################################################
  if(robust==FALSE)
    pcafun=wpca
  else
    pcafun=rwpca
  ##################################################################################
	for (i in 1:ep.n) 
  {
     if (DM.given)
        dist.vi<-dMat[,i]
     else
     {
        if (ep.given)
          dist.vi<-gw.dist(dp.locat, elocat, focus=i, p, theta, longlat) 
        else
          dist.vi<-gw.dist(dp.locat, focus=i, p=p, theta=theta, longlat=longlat) 
     }
		wt <-gw.weight(dist.vi,bw,kernel,adaptive)
		use <- wt>0
		wt <-wt[use]
    if(length(wt)<=5)
    {
       expr <- paste("Too small bandwidth at location: ", i)
       warning(paste(expr, "and the results can't be given there.", sep=", "))
       next
    }
		temp <- pcafun(x[use,],wt,nu=0,nv=k)
		w[i,,] <- temp$v
		d[i,] <- temp$d
    #Calculate the local scores
    if(scores)
    {
      scores.i <- c()
      for(j in 1:k)
      {
        score <- t(x[use,])*temp$v[,j]
        scores.i <- cbind(scores.i, apply(score,2,sum))
      }
      gwpca.scores[[i]] <- scores.i
    }
    
  }
	if (!is.null(rownames(x))) dimnames(w)[[1]] <- rownames(x)
	if (!is.null(colnames(x))) dimnames(w)[[2]] <- colnames(x)
	dimnames(w)[[3]] <- paste("PC",1:k,sep='')
	CV<-numeric(dp.n)
  if(cv)
     CV<-gwpca.cv.contrib(x,dp.locat,bw,k,robust,kernel,adaptive, p, theta, longlat,dMat)
	GW.arguments<-list(vars=vars,k=k, bw=bw, kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat, dp.n=dp.n, DM.given=DM.given,scores=scores)
	# And add and change this bit please #################################################
  if(robust==FALSE)
    d1 <- (d/(sum(wt)^0.5))^2
  else
    d1 <- d^2                     
  local.PV <- d1[, 1:k]/rowSums(d1) * 100
  var.names <- c()
  for(i in 1:k)
     var.names <- c(var.names, paste(paste("Comp", i, sep="."), "PV", sep="_"))
  win.var.pc1 <- max.col(abs(w[,,1]))
  res.df <- data.frame(local.PV, rowSums(local.PV), vars[win.var.pc1])
  names(res.df) <- c(var.names, "local_CP", "win_var_PC1")
  if (!is.null(polygons))
  {
    rownames(res.df) <- sapply(slot(polygons, "polygons"),
                                  function(i) slot(i, "ID"))
    SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=res.df,match.ID=F)
  }
  else
  {
    SDF <- SpatialPointsDataFrame(coords=elocat, data=res.df, proj4string=CRS(p4s), match.ID=F)
  }
  ################################################
  timings[["stop"]] <- Sys.time()
  res <- list(pca=pca.res,loadings = w, SDF=SDF,gwpca.scores=gwpca.scores, var=d1,  local.PV=local.PV, GW.arguments = GW.arguments, CV = CV, timings=timings)
  class(res) <-"gwpca"
  invisible(res) 
  ##################################################################################
 }
############################Layout function for outputing the GWPCA results
##Author: BL
print.gwpca<-function(x, ...)
{
  if(class(x) != "gwpca") stop("It's not a gwpca object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  vars <- x$GW.arguments$vars
  cat("\n   Variables concerned: ",vars)
  cat("\n   The number of retained components: ",x$GW.arguments$k)
  dp.n<- dim(x$loadings)[1]
  cat("\n   Number of data points:",dp.n)
  ################################################################ Print Linear
  cat("\n   ***********************************************************************\n")
  cat("   *                Results of Principal Components Analysis               *\n")
  cat("   ***********************************************************************\n")
  print(summary(x$pca,  loadings = TRUE, cutoff=0))
  
  #########################################################################
  cat("\n   ***********************************************************************\n")
    cat("   *   Results of Geographically Weighted Principal Components Analysis  *\n")
  cat("   ***********************************************************************\n")
  cat("\n   *********************Model calibration information*********************\n")
  cat("   Kernel function for geographically weighting:", x$GW.arguments$kernel, "\n")
  if(x$GW.arguments$adaptive)
    cat("   Adaptive bandwidth for geographically and temporally  weighting: ", x$GW.arguments$bw, " (number of nearest neighbours)\n", sep="")
  else
    cat("   Fixed bandwidth for geographically and temporally weighting: ", x$GW.arguments$bw, "\n")
  if (x$GW.arguments$DM.given)
    cat("   Distance metric for geographically weighting: A distance matrix is specified for this model calibration.\n")
  else
  {
    if (x$GW.arguments$longlat)
      cat("   Distance metric for geographically weighting: Great Circle distance metric is used.\n")
    else if (x$GW.arguments$p==2)
      cat("   Distance metric for geographically weighting: Euclidean distance metric is used.\n")
    else if (x$GW.arguments$p==1)
      cat("   Distance metric for geographically weighting: Manhattan distance metric is used.\n") 
    else if (is.infinite(x$GW.arguments$p))
      cat("   Distance metric for geographically weighting: Chebyshev distance metric is used.\n")
    else 
      cat("   Distance metric for geographically weighting: A generalized Minkowski distance metric is used with p=",x$GW.arguments$p,".\n")
    if (x$GW.arguments$theta!=0&&x$GW.arguments$p!=2&&!x$GW.arguments$longlat)
      cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GW.arguments$theta, "in radian.\n")   
  } 
  
  cat("\n   ****************     Summary of GWPCA information:    *****************\n")       
  var.names <- c()
  for(i in 1:x$GW.arguments$k)
     var.names <- c(var.names, paste("Comp", i, sep="."))
  cat("   Local variance: \n")
  local.SD <- t(apply(x$var[,1:x$GW.arguments$k], 2, summary))[,c(1:3,5,6)]
  rownames(local.SD) <- var.names
  if(x$GW.arguments$k==1) 
  { 
    local.SD <- matrix(local.SD, nrow=1)
    colnames(local.SD) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(local.SD) <- var.names
  }
  rnames<-rownames(local.SD)
  for (i in 1:length(rnames))
    rnames[i]<-paste("   ",rnames[i],sep="")
  rownames(local.SD) <-rnames 
  printCoefmat(local.SD)
  cat("   Local Proportion of Variance: \n")
  local.PV <-  t(apply(as(x$SDF, "data.frame")[,1:(x$GW.arguments$k+1), drop=FALSE], 2, summary))[,c(1:3,5,6)]
  if(x$GW.arguments$k==1) 
  { 
    local.PV <- matrix(local.SD, nrow=1)
    colnames(local.PV) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
  }
  rownames(local.PV) <-c(rnames, paste("   ","Cumulative",sep=""))
  printCoefmat(local.PV)
  cat("\n   ***********************************************************************\n")
  cat("   Program stops at:", as.character(x$timings$stop), "\n")
  invisible(x)
}
#######################################Bandwidth selection
# Function to find a fixed or adaptive bandwidth 
bw.gwpca<-function(data,vars,k=2, robust = FALSE,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F,dMat)
{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
  }
  else if (is(data, "data.frame")&&(!missing(dMat)))
     data<-data
  else
     stop("Given data must be a Spatial*DataFrame or data.frame object")
  data <- as(data, "data.frame")
  dp.n<-nrow(data)
  if (missing(dMat))
  {
    DM.given<-F
    if(dp.n <= 5000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
       stop("Dimensions of dMat are not correct")
  }
  if(missing(vars))
    stop("Variables input error")
  ##########Extract data
  col.nm<-colnames(data)
  var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if(length(var.idx)==0) stop("Variables input doesn't match with data")
  x<-data[,var.idx]
  x<-as.matrix(x)
  var.nms<-colnames(x)
  var.n<-ncol(x)
  #########Find the range of the fixed bandwidth
  if(adaptive)
  {
    upper<-dp.n
    lower<-2
  }
  else
  {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }
    ###!!!!!!!! Important note: if the distance matrix is not specified, the running time will be consuming very much by choosing the range of fixed bandwidth when the p is not 2; 
    ### because the range can't be decided by boundary box
    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
     
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
  # Add this bit please ############################################################
#  if(robust==FALSE)
#    pcafun=wpca
#  else
#    pcafun=rwpca
  ##################################################################################
    bw<-NA
    bw <- gold(gwpca.cv,lower,upper,adapt.bw=adaptive,x,dp.locat,k,robust,kernel,adaptive, p, theta, longlat,dMat) 
   # bw<-NA
#    if(approach=="cv"||approach=="CV")
#       bw <- optimize(bw.cv,lower=lower,upper=upper,maximum=FALSE,X=x,Y=y,kernel=kernel,
#       adaptive=adaptive, dp.locat=dp.locat, p=p, theta=theta, longlat=longlat,dMat=dMat,tol=.Machine$double.eps^0.25)
#    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
#       bw<-optimize(bw.aic,lower=lower,upper=upper,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
    bw

}
	
# Cross-validation function to optimally find a fixed or adaptive bandwidth... 
gwpca.cv <- function(bw,x,loc,k=2,robust=FALSE,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F,dMat) 
{
  if (missing(dMat))
    DM.given<-F
  else
    DM.given<-T
	n <- nrow(loc)
	m <- ncol(x)
	w <- array(0,c(n,m,k))
	#bw <- bw*bw
	score <- 0
  if(robust==FALSE)
    pcafun=wpca
  else
    pcafun=rwpca
	for (i in 1:n) 
  {
	   if (DM.given)
        dist.vi<-dMat[,i]
     else
     {
        dist.vi<-gw.dist(loc, focus=i, p=p, theta=theta, longlat=longlat) 
     }
		wt <-gw.weight(dist.vi,bw,kernel,adaptive)
		wt[i] <- 0
		use <- wt>0
		wt <-wt[use]
		if(length(wt)<=1)
    {
       score <-Inf
       expr <- paste("Too small bandwidth: ", bw)
       warning(paste(expr, "and the CV value can't be given there.", sep=", "))
       break
    }
		v <- pcafun(x[use,],wt,nu=0,nv=k)$v
		v <- v %*% t(v)
		score <- score + sum((x[i,] - x[i,] %*% v))**2
  }
  if(adaptive)
    cat("Adaptive bandwidth(number of nearest neighbours):", bw, "CV score:", score, "\n")
  else
    cat("Fixed bandwidth:", bw, "CV score:", score, "\n")
	score
}

# Contribution of each observation to the score statistic used in cross-validation for gwpca
# Outliers taken to correspond to high score (residual) values... 
gwpca.cv.contrib <- function(x,loc,bw, k=2,robust=FALSE,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F,dMat)
{
	if (missing(dMat))
    DM.given<-F
  else
    DM.given<-T
  n <- nrow(loc)
	m <- ncol(x)
	w <- array(0,c(n,m,k))
	score <- numeric(n)
  if(robust==FALSE)
    pcafun=wpca
  else
    pcafun=rwpca
	for (i in 1:n) 
  {
		if (DM.given)
        dist.vi<-dMat[,i]
     else
     {
        dist.vi<-gw.dist(loc, focus=i, p=p, theta=theta, longlat=longlat) 
     }
		wt <-gw.weight(dist.vi,bw,kernel,adaptive)
		wt[i] <- 0
		use <- wt>0
		wt <-wt[use]
		if(length(wt)<=1)
    {
       score[i] <-Inf
       expr <- paste("Too small bandwidth: ", bw)
       warning(paste(expr, "and the CV value can't be given there.", sep=", "))
       break
    }
		v <- pcafun(x[use,],wt,nu=0,nv=k)$v
		v <- v %*% t(v)
		score[i] <- sum((x[i,] - x[i,] %*% v))**2
  }
	score
 }
