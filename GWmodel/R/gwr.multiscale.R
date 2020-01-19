# Renamed as gwr.multiscale, i.e. multiscale GWR
# Geographically Weighted Regression with Parameter-Specific Distance Metrics  
# With Flexiable bandwidth selection
# Method (Yang, p44, 2014): select an optimum bandwidth for each coefficient(intercept inclusive),
# within each step of back-fitting, employing an existing criteria for a basic GWR

# Improvements: Set the corresponding threshold for each selected bandwidth, i.e. fix the bandwidth if the change of the selected bandwidths is less than threshold
gwr.multiscale <- function(formula, data, kernel="bisquare", adaptive=FALSE, criterion="dCVR", max.iterations=2000,threshold=0.00001, dMats, p.vals, theta.vals,longlat=FALSE,
                           bws0, bw.seled=rep(F, length(bws0)), approach = "AIC", bws.thresholds=rep(0.1, length(dMats)), bws.reOpts=5, verbose=F, hatmatrix=T, 
                           predictor.centered=rep(T, length(bws0)-1),nlower = 10)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    regression.points <- data
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
  #########Distance matrix is given or not
  dp.n <- nrow(dp.locat)
  ######Extract the data frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n <- ncol(x)
  #x2 <- x
  idx1 <- match("(Intercept)", colnames(x))
  ##################################################
  #####Linear regression
  lms <- lm(formula,data=data)
  #lms <-fastLm(formula,data=data) 
  lms$x <- x
  lms$y <- y
  #colnames(x)[1]<-"Intercept"
  #################################################
  ##########Centering and scaling predictors
  if(!is.na(idx1))
  {
    if(length(predictor.centered)!=var.n-1)
    {
      predictor.centered <- rep(T, length.out=var.n-1)
      warning("All the predictors will be centered, please check the parameter predictor.centered")
    }
    # if(length(predictor.scaled)!=var.n-1)
    # {
    #   predictor.scaled <- rep(F, length.out=var.n-1)
    #   warning("All the predictors will not be scaled, please check the parameter predictor.scaled")
    # }
  }
  else
  {
    if(length(predictor.centered)!=var.n)
    {
      predictor.centered <- rep(T, length.out=var.n)
      warning("All the predictors will be centered, please check the parameter predictor.centered")
    }
    # if(length(predictor.scaled)!=var.n-1)
    # {
    #   predictor.scaled <- rep(F, length.out=var.n-1)
    #   warning("All the predictors will not be scaled, please check the parameter predictor.scaled")
    # }
  }
  
  n.cent <- length(which(predictor.centered))
  if(!is.na(idx1))
  {
    colnames(x)[idx1]<-"Intercept"
    x1 <- x[,-idx1]
  }
  else
  {
    x1 <- x
  }
  #n.scal <- length(which(predictor.centered))
  if(n.cent>1)
    predictors.centered.means <- colMeans(x1[,predictor.centered])
  else if (n.cent==1)
    predictors.centered.means <- mean(x1[,predictor.centered])
  else
    predictors.centered.means <- NA
  if(n.cent>=1)
  {
    #meancent <- partial(scale, scale=FALSE)
    x1[,predictor.centered] <- scale(x1[,predictor.centered], scale=FALSE)
    
  }
  if(!is.na(idx1))
  {
    x1 <- cbind(1, x1)
  }
  colnames(x1) <-  colnames(x)
  
  #####Centering predictors
  #####Centered predictors x1
  #################################################
  
  ##Calculate the initial bandwidths: find an optimum bandwidth for each univariate model (dependent variable~each independent variable)
  allvars <- all.vars(formula)
  DeVar <- allvars[1]
  InDevars <- colnames(x)
  ##########################Check the distance matrices
  ###check the distance matrices and bandwidths
  if(missing(dMats))
  {
    dMats <- list()
    if(missing(p.vals))
    {
      p.vals <- 2
      dMat <- gw.dist(dp.locat,longlat=longlat)
      for(i in 1:var.n)
        dMats[[i]] <- dMat
    }
    else
    {
      p.vals <- rep_len(p.vals, var.n)
      if(missing(theta.vals))
        theta.vals <- rep_len(0, var.n)
      else
        theta.vals<- rep_len(theta.vals, var.n)
      for(i in 1:var.n)
      {
        dMats[[i]] <- gw.dist(dp.locat, p=p.vals[i], theta=theta.vals[i],longlat=longlat)
      }
    }
  }
  else if(is.list(dMats))
  {
    if(length(dMats)!=var.n)
      stop("Please specify a distance matrix for each independent variable!")
    else
    {
      for(i in 1:var.n)
      {
        dim.dMat<-dim(dMats[[i]])
        if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
          stop("Dimensions of dMat are not correct")
      }
    }
  }
  else
  {
    stop("dMats are not correct")
  }  
  ############################### Intialize the bandwidth
  if(missing(bws0))
  {
    bws0 <- numeric(var.n)
    cat("------ Calculate the initial bandwidths for each independent variable ------\n")
    for(i in 1:var.n)
    {
      if(InDevars[i]=="Intercept")
        fml <- Generate.formula(DeVar,c("1"))
      else
      {
        fml <- Generate.formula(DeVar,c(InDevars[i]))
        #fml <- paste(fml, "1", sep="-")
      }
      cat("Now select an optimum bandwidth for the model: ", fml,"\n")
      part1<-paste("bws0[i]<-bw.gwr(",fml,sep="")
      part2<-"data=regression.points,kernel=kernel,approach=approach,dMat=dMats[[i]])"
      expression<-paste(part1,part2,sep=",")
      print(expression)
      eval(parse(text=expression))
    }
    cat("------            The end for the initial selections              ------\n")
  }
  else
  {
    bws0 <- rep(bws0, length.out=var.n)
  }
  #################
  if(length(bw.seled)!=var.n)
    bw.seled <- rep(F, length.out=var.n)

  cat("------   Calculate the initial beta0 from the above bandwidths    ------\n")
  #dMat <- gw.dist(dp.locat=dp.locat, p=2, theta=0, longlat=longlat)
  dMat <- dMats[[1]]
  bw.int0 <- bw.gwr1(x1, y, dp.locat,approach=approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
  #betas <- gwr.q(x, y, dp.locat, adaptive=adaptive, bw=bw.int0, kernel=kernel,dMat=dMat)
  #######Re-initialize the betas by a simple back-fitting process
  #betas <- gwr.backfit(x, y, betas,dp.locat,dp.locat, FALSE,criterion, adaptive, bws0, kernel,dMats, max.iterations,threshold*100)
  #####Hatmatrix for the whole process
  if(hatmatrix)
  {
    Shat <- matrix(nrow=dp.n, ncol=dp.n)
    S.arrays <- array(dim=c(var.n, dp.n, dp.n))
    C <- array(dim=c(dp.n,var.n,dp.n))
  }
  
  res <- gwr.q2(x1, y, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bw.int0, kernel=kernel,dMat=dMat)
  betas <- res[[1]]
  if(hatmatrix)
  {
    Shat <- res[[2]]
    C <- res[[3]]
    idm <- diag(var.n)
    for(i in 1:var.n)
      for(j in 1:dp.n)
      {
        S.arrays[i,j,] <- x1[j,i]*(idm[i,]%*%C[j,,])
      }
  }
  #betas <- gwr.backfit(x, y, betas,dp.locat,dp.locat, FALSE,criterion, adaptive, bws0, kernel,dMats, max.iterations,0.0001)
  cat("------            The end for calculating the initial beta0              ------\n") 
  cat("------ Select the optimum bandwidths for each independent variable via ", approach, " aproach ------\n")
  ieration <-0 
  #bws1 <- bws0
  bws.vars <- bws0
  bws.change.NO <- numeric(var.n)
  criterion.val <- 10000000  
  #yhat.i <- betas*x
  #print(yhat.i)
  resid.i <- y - gw.fitted(x1, betas)
  RSS0 <- sum(resid.i^2)
  RSS1 <- 0
  RSS.vals <- c(RSS0, RSS1, criterion.val)
  cat("*****  The back-fitting process for model calibration with bandwiths selected *****\n")
  while((ieration < max.iterations) && criterion.val > threshold)
  { 
    #AICcs <- numeric(var.n)
    cat("    Iteration ", ieration+1, ":\n")   
    for(i in 1:var.n)
    {
      dMat <- dMats[[i]]  
      f.i <- betas[,i]*x1[,i]
      y.i <- resid.i + f.i
      if(bw.seled[i])
        bw.i <- bws0[i]
      else
      {
        cat("Now select an optimum bandwidth for the variable: ", InDevars[i], "\n")
        bw.i <- bw.gwr1(matrix(x1[, i], ncol = 1), y.i, dp.locat, approach = approach, kernel = kernel, adaptive = adaptive, dMat, verbose = verbose, nlower = nlower)
        cat("The newly selected bandwidth for variable ", InDevars[i])
        cat(" is: ", bw.i, "\n")
        cat("The bandwidth used in the last ieration is: ", bws0[i])
        cat(" and the difference between these two bandwidths is: ", abs(bw.i - bws0[i]), "\n")
        if (abs(bw.i - bws0[i]) > bws.thresholds[i]) {
          cat("The bandwidth for variable ", InDevars[i])
          cat(" will be continually selected in the next ieration.\n")
          bws.change.NO[i] <- 0
        }
        else {
          bws.change.NO[i] <- bws.change.NO[i] + 1
          if(bws.change.NO[i] < bws.reOpts)
          {
            cat("The bandwidth for variable ", InDevars[i])
            cat(" seems to be converged for ", bws.change.NO[i], " times.")
            cat("It will be continually optimized in the next ", bws.reOpts-bws.change.NO[i], " times\n")
          }
          else
          {
            cat("The bandwidth for variable ", InDevars[i])
            cat(" seems to be converged and will be kept the same in the following ierations.\n")
            bw.seled[i] <- T
          }
        }
      }
      bws0[i] <- bw.i       
      res <- gwr.q2(matrix(x1[,i], ncol=1), y.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bw.i, kernel=kernel,dMat=dMat)
      betai <- res[[1]]
      if(hatmatrix)
      {
        Si <- res[[2]]
        ###See Yu et al. 2018, eq. 18 on P8
        S.arrayi <- S.arrays[i,,]
        S.arrays[i,,] <- Si%*%S.arrayi + Si - Si%*%Shat
        Shat <- Shat- S.arrayi + S.arrays[i,,] 
      }
      
      #betai <- gwr.q(matrix(x[,i], ncol=1), y.i, loc=dp.locat, adaptive=adaptive, bw=bw.i, kernel=kernel,dMat=dMat)
      #AICcs[i] <- gwr.aic(bw.i, matrix(x[,i], ncol=1), y.i, kernel, adaptive, dp.locat, dMat=dMat, verbose=F)
      #yhat.i[,i] <- betai*x[,i]
      betas[,i] <- betai
      resid.i <- y - gw.fitted(x1, betas)
      #resid.i <- ehat(y.i, matrix(x[,i], ncol=1), matrix(betas[,i], ncol=1))
      #betas[,i] <- betai
    }
    bws.vars <- rbind(bws.vars, bws0)
    #AICc.vals <- rbind(AICc.vals, AICcs) 
    RSS1 <- sum((y - gw.fitted(x1, betas))^2)   
    if(criterion=="CVR")
    {
      criterion.val <- abs(RSS1-RSS0)
      cat("    Ieration ", ieration, "the change value of RSS (CVR) is: ", criterion.val,"\n")
    }
    else
    {
      criterion.val <- sqrt(abs(RSS1-RSS0)/RSS1)
      cat("    Ieration ", ieration+1, "the differential change value of RSS (dCVR) is: ", criterion.val,"\n") 
    }
    RSS0 <- RSS1  
    cat("----------End of    Iteration ", ieration+1, "----------\n") 
    RSS.vals <- rbind(RSS.vals, c(RSS0, RSS1, criterion.val)) 
    ieration <- ieration+1            
  }
  #####Output
  yhat <- gw.fitted(x1, betas)
  residual <- y - yhat
  ################
  RSS.gw <- RSS1
  sigma.hat21 <- RSS.gw/dp.n
  yss.g <- sum((y - mean(y))^2)
  R2.val <-  1-RSS.gw/yss.g
  R2adj <- NA
  AICc <- NA
  if(hatmatrix)
  {
    tr.Shat <- sum(diag(Shat))
    tr.StShat<-sum(Shat^2)
    edf<- dp.n - 2*tr.Shat + tr.StShat
    AICc <- dp.n*log(sigma.hat21) + dp.n*log(2*pi) + dp.n *((dp.n + tr.Shat) / (dp.n - 2 - tr.Shat))
    R2adj <- 1-(1-R2.val)*(dp.n-1)/(edf-1)
  }
  
  #sigma.hat11<-RSS.gw/(dp.n-2*tr.Shat+tr.StShat)
  GW.diagnostic<-list(RSS.gw=RSS.gw,AICc=AICc,R2.val=R2.val, R2adj = R2adj)
  #########################Retrive the intercept
  if(!is.na(idx1) && n.cent>=1)
  {
    idx.cent <- which(predictor.centered)
    betas.cent <- c()
    for(i in 1:n.cent)
      betas.cent <- cbind(betas.cent, betas[,idx.cent[i]+1]*predictors.centered.means[i])
    beta0 <- betas[,idx1] - apply(betas.cent,1, sum)
    betas[,idx1] <- beta0
  }
  ############################
  vdgwr.df <- data.frame(betas, yhat, residual)
  #colnames(vdgwr.df) <- c(colnames(x), "yhat", "residual",paste(colnames(betas), "SE", sep="_"), paste(colnames(betas), "TV", sep="_"))
  colnames(vdgwr.df) <- c(colnames(x), "yhat", "residual")
  griddedObj <- F
  if (is(regression.points, "Spatial"))
  { 
    if (is(regression.points, "SpatialPolygonsDataFrame"))
    {
      polygons<-polygons(regression.points)
      #SpatialPolygons(regression.points)
      #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
      #  function(i) slot(i, "ID"))
      SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=vdgwr.df,match.ID=F)
    }
    else
    {
      griddedObj <- gridded(regression.points)
      SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwr.df, proj4string=CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj 
    }
  }
  else
    SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwr.df, proj4string=CRS(p4s), match.ID=F)  
  #GW.arguments<-list(formula=formula,rp.given=T,hatmatrix=T,criterion="RSS",bws=bws0, kernel=kernel,adaptive=adaptive)  
  timings[["stop"]] <- Sys.time()
  GW.arguments<-list(formula=formula,criterion=criterion,bws=bws0, kernel=kernel,adaptive=adaptive, hatmatrix=hatmatrix)  
  #res <- list(SDF=SDF, GW.arguments=GW.arguments,lm=lms,timings=timings,this.call=this.call)
  #class(res) <- "psdmgwr"
  #invisible(res)
  res <- list(SDF=SDF,GW.arguments=GW.arguments,GW.diagnostic=GW.diagnostic,lm=lms, bws.vars,timings=timings,this.call=this.call)
  class(res) <- "multiscalegwr"
  invisible(res) 
}

############################Layout function for outputing the Multiscale(PSDM) GWR results
##Author: BL
print.multiscalegwr<-function(x, ...)
{
  if(class(x) != "multiscalegwr") stop("It's not a multi-scale gwr object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GW.arguments$formula)
  var.n<-length(x$lm$coefficients)
	cat("\n   Dependent (y) variable: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	dp.n<-length(x$lm$residuals)
	cat("\n   Number of data points:",dp.n)
	#########################################################################
	cat("\n   ***********************************************************************\n")
    cat("   *             GWR with Parameter-Specific Distance Metrics            *\n")
	cat("   ***********************************************************************\n")
	cat("\n   *********************Model calibration information*********************\n")
	cat("   Kernel function:", x$GW.arguments$kernel, "\n")
	if(x$GW.arguments$adaptive)
	   cat("   Adaptive bandwidths for each coefficient(number of nearest neighbours): \n") 
  else
     cat("   Fixed bandwidths for each coefficient: \n")
	bws <- matrix(x$GW.arguments$bws,nrow=1)
  rownames(bws) <- c("   Bandwidth ")
  colnames(bws) <- names(x$lm$coefficients)
  printCoefmat(bws)
	cat("\n   ****************Summary of GWR coefficient estimates:******************\n")       
		df0 <- as(x$SDF, "data.frame")[,1:var.n, drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	if(var.n==1) 
    { 
      CM <- matrix(CM, nrow=1)
      colnames(CM) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
      rownames(CM) <- names(x$SDF)[1]
    }
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM)
	cat("   ************************Diagnostic information*************************\n")
	#rownames(diag.mat) <- c("   Residual sum of squares", "   AICc value:","   R-square value","   Adjusted R-square value")
	options("scipen"=999)
	cat("   Residual sum of squares: ", x$GW.diagnostic$RSS.gw, "\n")
	cat("   R-square value: ", x$GW.diagnostic$R2.val, "\n")
	if(x$GW.arguments$hatmatrix)
	{
	  cat("   Adjusted R-square value: ", x$GW.diagnostic$R2adj, "\n")
	  cat("   AICc value: ", x$GW.diagnostic$AICc, "\n")
	}
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}

##Back-fitting algorithm, see Wenbai's thesis
gwr.backfit <- function(x, y, betas,dp.locat,rp.locat,hatmatrix, criterion="CVR", adaptive=F, bws, kernel="bisquare",dMats, max.iterations=100,threshold=0.0001)
{
   ieration <-0
   y.i <- y
   var.n <- ncol(x)
   dp.n <- nrow(x)
   rp.n <- nrow(rp.locat)
   criterion.val <- 10000000
   resid.i <- y - gw.fitted(x, betas)
   #resid.i <- ehat(y, x, betas)
   RSS0 <- sum(resid.i^2)
   RSS1 <- 0
   #if(criterion=="RSS")
   #{
     #
     #
   
   #  RSS0 <- rss(y, x, betas)
  #   
   #}
   #else
   #{
    #  betas0 <- betas  
   #}  
   cat("*************The back-fitting process*************\n")
   while((ieration < max.iterations) && criterion.val > threshold)
   {    
     for(i in 1:var.n)
     {
       bw <- bws[i]
       dMat <- dMats[[i]] 
       f.i <- betas[,i]*x[,i]
       y.i <- f.i+resid.i
       betai <- gwr.q(matrix(x[,i], ncol=1), y.i, loc=dp.locat, adaptive=adaptive, bw=bw, kernel=kernel,dMat=dMat)
       
       resid.i <- y.i - betai*x[,i]
       #resid.i <- ehat(y.i, matrix(x[,i], ncol=1), matrix(betas[,i], ncol=1))
       betas[,i] <- betai
     }
     #RSS1 <- rss(y,x,betas)
     RSS1 <- sum((y - gw.fitted(x, betas))^2)   
     if(criterion=="CVR")
     {
         #RSS1 <- sum((y - gwr.fitted(x, betas))^2)
         #cat("RSS1: ", RSS1,"\n")
         #cat("RSS0: ", RSS0,"\n")
         criterion.val <- abs(RSS1-RSS0)
         cat("    Ieration ", ieration, "the change value of RSS (CVR) is: ", criterion.val,"\n")
     }
     else
     {
          criterion.val <- sqrt(abs(RSS1-RSS0)/RSS1)
          cat("    Ieration ", ieration, "the differential change value of RSS (dCVR) is: ", criterion.val,"\n") 
     }
     RSS0 <- RSS1  
     ieration <- ieration+1            
   }
   #calculate the hatmatrix
   betas
}
###For bandwidth selection
bw.gwr2<-function(x, y, dp.locat,approach="AIC",kernel="gaussian",adaptive=FALSE,dMat, verbose=F, nlower = 10)
{
  dp.n <-  nrow(dp.locat)
  if(adaptive)
  {
    upper <- dp.n
    lower <- nlower
  }
  else
  {
      upper<-range(dMat)[2]
      lower<-upper/5000
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA
    if(approach=="cv"||approach=="CV")
       bw <- gold(gwr.cv,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat,dMat=dMat,verbose=verbose)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw <- gold(gwr.aic1,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat,dMat=dMat,verbose=verbose)
    bw
}
gwr.cv2 <- function(bw, X, Y, kernel,adaptive, dp.locat,dMat)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   betas <- matrix(nrow=dp.n,ncol=var.n) 
  CV<-numeric(dp.n)
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    betas[i,] <- try(fun1(X,Y,W.i))
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))
    W.i[i]<-0
    gw.resi<-try(fun1(X,Y,W.i))
    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  CV.score<-t(CV) %*% CV
  res<- list(CV.score, betas)
}
gwr.q2 <- function(x, y, loc, adaptive=F, hatmatrix,bw=sqrt(var(loc[,1])+var(loc[,2])),
                  kernel, p, theta, longlat,dMat, wt2=rep(1,nrow(loc)))
{
  if (missing(dMat))
     DM.given <- F
  else
     DM.given <- T
  dp.n <- nrow(loc)
  var.n <- ncol(x)
  betas <- matrix(nrow=dp.n, ncol=var.n)
  S <- matrix(nrow=dp.n, ncol=dp.n)
  C <- array(dim=c(dp.n, var.n, dp.n))
  for (i in 1:dp.n)
  {
    if(DM.given)
       dist.vi <- dMat[,i]
    else
       dist.vi <- gw.dist(loc, focus=i, p, theta, longlat)
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw_reg(x,y,as.vector(W.i*wt2),hatmatrix=hatmatrix,i)
    betas[i,]<-gw.resi[[1]]
    if(hatmatrix)
    {
      S[i,]<-gw.resi[[2]]
      C[i,,]<-gw.resi[[3]]
    }
  }
  colnames(betas) <- colnames(x)
  res <- list(betas, S,C)
  res
}

####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
gwr.aic1<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   #########Distance matrix is given or not

  if (is.null(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  S<-matrix(nrow=dp.n,ncol=dp.n)
  betas <-matrix(nrow=dp.n, ncol=var.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    res<- try(gw_reg(X,Y,W.i,TRUE,i))
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    #fun2<-function(X,W.i) {Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}}
    #Ci<-try(fun2(X,W.i))
    
    #Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
   # gw.resi<-gw.reg(X,Y,W.i,hatmatrix=T,focus=i)
    #betas[i,]<-gw.resi[[1]] ######See function by IG
    #S[i,]<-gw.resi[[2]]
    if(!inherits(res, "try-error"))
    {
      S[i,]<-res[[2]]   
      betas[i,] <- res[[1]]
    }
    else
    {
      S[i,]<-Inf
      break
    }  
  }
  
  if (!any(is.infinite(S)))
  {
    #tr.S<-sum(diag(S))
#    RSS.gw<-t(Y)%*%t(diag(dp.n)-S)%*%(diag(dp.n)-S)%*%Y
#    sigma.hat2 <- RSS.gw/dp.n
#    AICc<-dp.n*log(sigma.hat2) + dp.n*log(2*pi) + dp.n *((dp.n + tr.S) / (dp.n - 2 - tr.S))
     AICc<-AICc(Y,X,betas, S)
  }
  else
    AICc<-Inf
  if(verbose)
  {     
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
    else
      cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  }
  if(is.nan(AICc))
      AICc <- Inf
  AICc
}




















