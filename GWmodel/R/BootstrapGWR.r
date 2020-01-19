#Bootstrap GWR code edited from Harry's code
#Basic test statistic, modified test statistic and localised test statistic are returned
#R: number of random samples
gwr.bootstrap <- function(formula, data, kernel="bisquare",approach="AIC", R=99,k.nearneigh=4,adaptive=FALSE, p=2, theta=0, longlat=FALSE,dMat,verbose=FALSE)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  n.sim.rep <- 10
  ######data points
  ###Create the adjency matrix for calculating the ERR, SMA and LAG model
  polygons <- NULL
  if (is(data, "SpatialPolygonsDataFrame"))
  {
     polygons <- polygons(data)
	 dp.locat <- coordinates(data)
	 gnb <- poly2nb(polygons)
	 glw <- nb2listw(gnb)
	 W.adj <- listw2mat(glw) 
  }    
  else if(is(data, "SpatialPointsDataFrame"))
  {
    dp.locat <- coordinates(data)
	gnb <- knn2nb(knearneigh(dp.locat, k=k.nearneigh), sym = T)
	glw <- nb2listw(gnb)
	W.adj <- listw2mat(glw) 
	#griddedObj <- gridded(dp.locat)
  }
  else
    stop("Given regression data must be Spatial*DataFrame")
  sp.data <- data
  data <- as(data, "data.frame")
  ####################################################GWR
	  #########Distance matrix is given or not
  dp.n <- nrow(dp.locat)
  if (missing(dMat))
  {
      dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
  }
  else
  {
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
       stop("Dimensions of dMat are not correct")
  } 
  ###################
	#####LM model
    ols.model <- lm(formula, data)
    err.model <- errorsarlm(formula,data,listw=glw,method='spam')
	dep.var <- deparse(eval(err.model$call$formula)[[2]])
    sma.model <- spautolm(formula,data,listw=glw,family='SMA')
    sma.model$model <- ols.model$model # As sma returns NULL otherwise
    lag.model <- lagsarlm(formula,data,listw=glw,method='spam')
  
     ###Basic GWR model
   bw <- bw.gwr3(formula,data=sp.data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	 gwr.model <- gwr.basic(formula,data=sp.data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat) 
   #####
	# For modified test statistic
    ols.bst <- parametric.bs(ols.model,dep.var,dp.locat,W.adj,gwrtvar,R=R, report=n.sim.rep,formula=formula, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	err.bst <- parametric.bs(err.model,dep.var,dp.locat,W.adj,gwrtvar,R=R, report=n.sim.rep,formula=formula, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)	
	sma.bst <- parametric.bs(sma.model,dep.var,dp.locat,W.adj,gwrtvar,R=R, report=n.sim.rep,formula=formula, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
    lag.bst <- parametric.bs(lag.model,dep.var,dp.locat,W.adj,gwrtvar,R=R, report=n.sim.rep,formula=formula, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
    sp.data <- SpatialPointsDataFrame(dp.locat,data, match.ID=FALSE)
	actual.t = gwrtvar(sp.data, formula,approach, kernel, adaptive,dMat,verbose=verbose) 
	results.t = rbind(ci.bs(ols.bst,0.95),pval.bs(ols.bst,actual.t),				  
				  ci.bs(err.bst,0.95),pval.bs(err.bst,actual.t),
				  ci.bs(sma.bst,0.95),pval.bs(sma.bst,actual.t),
                  ci.bs(lag.bst,0.95),pval.bs(lag.bst,actual.t))
	rownames(results.t) = c("   Modified statistic for MLR at 95% level","   p value to accept null hypothese(MLR)","   Modified statistic for ERR at 95%","   p value to accept null hypothese (ERR)",
	"   Modified statistic for SMA at 95% level","   p value to accept null hypothese (SMA)","   Modified statistic for LAG at 95% level","   p value to accept null hypothese (LAG)")
	###Localized test statistic
	err.bsm <- parametric.bs.local(err.model,dep.var,dp.locat,W.adj,gwrt.err,R=R,report=n.sim.rep,formula=formula, glw, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	mlr.bsm <- parametric.bs.local(ols.model,dep.var,dp.locat,W.adj,gwrt.mlr,R=R,report=n.sim.rep,formula=formula, approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	sma.bsm <- parametric.bs.local(sma.model,dep.var,dp.locat,W.adj,gwrt.sma,R=R,report=n.sim.rep,formula=formula, glw,approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
  lag.bsm <- parametric.bs.local(lag.model,dep.var,dp.locat,W.adj,gwrt.lag,R=R,report=n.sim.rep,formula=formula, glw,approach=approach, kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
  actual.m.err <- gwrt.err(sp.data,formula,glw,approach, kernel, adaptive,dMat,verbose=verbose)
	actual.m.mlr <- gwrt.mlr(sp.data,formula,approach, kernel, adaptive,dMat,verbose=verbose)
	actual.m.sma <- gwrt.sma(sp.data,formula,glw,approach, kernel, adaptive,dMat,verbose=verbose)
  actual.m.lag <- gwrt.lag(sp.data,formula,glw,approach, kernel, adaptive,dMat,verbose=verbose)
	indep.vars <- names(ols.model$coefficients)
	var.n <- length(indep.vars)
	idx1 <- match("(Intercept)",indep.vars)
    if(!is.na(idx1))
      indep.vars[idx1] <-"Intercept"
    colnames(results.t) <- indep.vars	  
	err.t.local <- c()
	mlr.t.local <- c()
	sma.t.local <- c()
  lag.t.local <- c()
	err.p.local <- c()
	mlr.p.local <- c()
	sma.p.local <- c()
  lag.p.local <- c()
	for(i in 1:var.n)
	{
	  err.q.vec <- c()
	  mlr.q.vec <- c()
	  sma.q.vec <- c()
    lag.q.vec <- c()
	  for(j in 1:R)
	  {
	    err.q.vec <- rbind(err.q.vec, err.bsm[[j]][,i])
	    mlr.q.vec <- rbind(mlr.q.vec, mlr.bsm[[j]][,i])
	    sma.q.vec <- rbind(sma.q.vec, sma.bsm[[j]][,i])
      lag.q.vec <- rbind(lag.q.vec, lag.bsm[[j]][,i])
	  }
	  err.t.local <- cbind(err.t.local, ci.bs(err.q.vec,0.95))
	  mlr.t.local <- cbind(mlr.t.local, ci.bs(mlr.q.vec,0.95))
	  sma.t.local <- cbind(sma.t.local, ci.bs(sma.q.vec,0.95))
    lag.t.local <- cbind(lag.t.local, ci.bs(lag.q.vec,0.95))
	  err.p.local <- cbind(err.p.local, pval.bs(err.q.vec,actual.m.err[,i]))
	  mlr.p.local <- cbind(mlr.p.local, pval.bs(mlr.q.vec,actual.m.mlr[,i]))
	  sma.p.local <- cbind(sma.p.local, pval.bs(sma.q.vec,actual.m.sma[,i]))
    lag.p.local <- cbind(lag.p.local, pval.bs(lag.q.vec,actual.m.lag[,i]))
	}
	

  #print()
  
  bsm.local.df <- data.frame((gwr.model$SDF@data)[,c(1:(var.n+2), (var.n+6):(var.n*2+5))], mlr.p.local, err.p.local, sma.p.local, lag.p.local)
  gwr.nms <- names(gwr.model$SDF@data)
	names(bsm.local.df) <- c(gwr.nms[c(1:(var.n+2), (var.n+6):(var.n*2+5))], paste(indep.vars, "MLR_p", sep="_"), paste(indep.vars, "ERR_p", sep="_"),
	                         paste(indep.vars, "SMA_p", sep="_"),paste(indep.vars, "LAG_p", sep="_"))	
	if (!is.null(polygons))
    {
      # polygons<-polygons(data)
       #SpatialPolygons(regression.points)
       rownames(bsm.local.df) <- sapply(slot(polygons, "polygons"),
                            function(i) slot(i, "ID"))
       SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=bsm.local.df,match.ID=F)
    }
    else
    {
       SDF <- SpatialPointsDataFrame(coords=dp.locat, data=bsm.local.df, proj4string=CRS(p4s), match.ID=F)
    }
    timings[["stop"]] <- Sys.time()
	res <- list(formula=formula, results = results.t, SDF=SDF,timings=timings,this.call=this.call)
	class(res) <-"gwrbsm"
    invisible(res)
}
######################
print.gwrbsm <- function(x, ...)
{
  if(class(x) != "gwrbsm") stop("It's not a gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$formula)
	cat("\n   Dependent (y) variable: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	dp.n<-nrow(x$SDF@data)
	cat("\n   Number of data points:",dp.n)
 indep.vars <- colnames(x$results)
 	  var.n <- length(indep.vars)
################################################################ Print Global statistic
 	cat("\n   ***********************************************************************\n")
	  cat("   *                             Bootstrap GWR                           *\n")
      cat("   ***********************************************************************\n")
	  cat("   ***                 Geographically weighted regression              ***\n")
    df0 <- x$SDF@data[, 1:(var.n*2+2)]
    if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	  CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM) 
	  cat("   ***********************************************************************\n")
	  cat("   ***                      Modified test statistic                    ***\n")
      cat("   *Comparison with a multiple linear regression model (MLR):\n\n")
	  
 	  cat("    Modified statistic for MLR at 95% level:\n")
      dm <- matrix(x$results[1,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n    p value to accept null hypothese(MLR):\n")
	  dm <- matrix(x$results[2,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n   *Comparison with a simultaneous autoregressive error model (ERR):\n\n")
 	  cat("    Modified statistic for ERR at 95%: \n")
	  dm <- matrix(x$results[3,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n    p value to accept null hypothese(ERR):\n")
	  dm <- matrix(x$results[4,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n   *Comparison with a moving average error model (SMA):\n\n")
 	  cat("    Modified statistic for SMA at 95%: \n")
	  dm <- matrix(x$results[5,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n    p value to accept null hypothese(SMA): \n")
	  dm <- matrix(x$results[6,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n   *Comparison with a simultaneous autoregressive lag model (LAG):\n\n")
 	  cat("    Modified statistic for LAG at 95%: \n")
	  dm <- matrix(x$results[7,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
	  cat("\n    p value to accept null hypothese(LAG): \n")
	  dm <- matrix(x$results[8,], nrow=1)
	  rownames(dm) <- "   "
	  colnames(dm) <- indep.vars
	  printCoefmat(dm)
    cat("   ***********************************************************************\n")
	  cat("   ***                      Localized test statistic                   ***\n")
	  df1 <- x$SDF@data[, (1:(var.n*2+2))*(-1)]
        if (any(is.na(df0))) {
            df1 <- na.omit(df1)
            warning("NAs in coefficients dropped")
        }
	  CM <- t(apply(df1, 2, summary))[,c(1:3,5,6)]
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM)
	cat("   *** Note the '_p' means the p value from the localised pseudo-t statistic.\n")
	cat("   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}
	
################################################################################

# Bootstrapping functions - source('bootpara_ig_ph.R')

generate.lm.data <- function(obj,W,dep.var) {
	generate.data.lm <- function(obj,W,dep.var) {
		x = obj$model
		#yname = deparse(obj$terms[[2]])
		yname <- dep.var
		ypos = match(yname,colnames(x))
		y = fitted(obj) + rnorm(nrow(x),0,summary(obj)$sigma)
		x[,ypos] = y
		data.frame(x)}
# Error SAR lm
	generate.data.errorsarlm <- function(obj,W,dep.var) {
		x = as.matrix(cbind(obj$y,obj$X[,-1]))
		# x = as.matrix(obj$lm.model$model)
		#yname <- deparse(eval(obj$call$formula)[[2]])
		yname <- dep.var
		ypos=1
		# ypos = match(yname,colnames(x))
		y = fitted(obj)+solve(diag(nrow(x))-obj$lambda*W,rnorm(nrow(x),0,sqrt(obj$s2)))
		x[,ypos] = y
		colnames(x)[ypos] <- yname
		x <- data.frame(x)
		x}
# Error SMA lm
	generate.data.smalm <- function(obj,W,dep.var) {
		x = as.matrix(obj$model)
		yname <- dep.var
		#yname <- all.vars(eval(obj$call$formula))[1]
		ypos = match(yname,colnames(x))
		y = fitted(obj)+(diag(nrow(x))+obj$lambda*W) %*% rnorm(nrow(x),0,sqrt(obj$fit$s2))
		x[,ypos] = y
		x <- data.frame(x)
		}
# Lagged SAR lm
	generate.data.lagsarlm <- function(obj,W,dep.var) {
		x = as.matrix(cbind(obj$y,obj$X[,-1]))
		yname <- dep.var
		#yname <- all.vars(eval(obj$call$formula))[1]
		ypos=1
		y = solve(diag(nrow(x))-obj$rho*W,obj$fitted+rnorm(nrow(x),0,sqrt(obj$s2)))
		x[,ypos] = y
		colnames(x)[ypos] = yname
		x <- data.frame(x)}
# What kind of model is it
	model.type <- function(obj) {
		if (class(obj) == "lm") return("lm")
		if (class(obj) == "spautolm") return("spautolm")
		if (class(obj) != "sarlm") stop("Unsupported regression type.")
		if (obj$type=="error") return("errorsarlm")
		if (obj$type=="lag") return("lagsarlm") }
# Do the simulation
	switch(model.type(obj),
		lm=generate.data.lm(obj,W,dep.var),
		lagsarlm=generate.data.lagsarlm(obj,W,dep.var),
		errorsarlm=generate.data.errorsarlm(obj,W,dep.var), 
		spautolm=generate.data.smalm(obj,W,dep.var))
		}

parametric.bs <- function(obj,dep.var,dp.locat,W,bsfun,R=100,report=NULL,...) {
	result = NULL
	for (i in 1:R) 
	{
		if (! is.null(report) & (i %% report == 0)) 
			cat(sprintf("Iteration %5d\n",i)) 
		dset <- generate.lm.data(obj,W,dep.var)
		sp.dset <- SpatialPointsDataFrame(dp.locat,dset, match.ID=FALSE)
		#print(class(sp.dset))
		result <- rbind(result,bsfun(sp.dset,...))
	}
	result }
####Localized statistic	
parametric.bs.local <- function(obj,dep.var,dp.locat,W,bsfun,R=100,report=NULL,...) {
	result <- list()
	for (i in 1:R) 
	{
		if (! is.null(report) & (i %% report == 0)) 
			cat(sprintf("Iteration %5d\n",i)) 
		dset <- generate.lm.data(obj,W,dep.var)
		sp.dset <- SpatialPointsDataFrame(dp.locat,dset, match.ID=FALSE) 
		result[[i]] <- bsfun(sp.dset,...)
	}
	result }
	
####Tiny functions
se.bs <- function(bs.out) apply(bs.out,2,sd)
bias.bs <- function(bs.out,stat) apply(bs.out,2,mean) - stat
ci.bs <- function(bs.out,ci) apply(bs.out,2,quantile,ci)
pval.bs <- function(bs.out,stat) apply(sweep(bs.out,2,stat,'>'),2,sum)/(nrow(bs.out)+1)

# Modified test statistic
gwrtvar <- function(data,formula, approach, kernel, adaptive,dMat,verbose=FALSE) {
	bw <- bw.gwr3(formula,data=data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	gwr.model <- gwr.basic(formula,data=data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat)
	var.n <- length(gwr.model$lm$coefficients)
	coefs <- gwr.model$SDF@data[,1:var.n] # Key
	sds   <- gwr.model$SDF@data[,(var.n+6):(var.n*2+5)] # Key and corrected from CB & IG
	tvals <- coefs/sds
	apply(tvals,2,sd)}
###########################localised test statistic
gwrt.mlr <- function(data,formula,approach, kernel, adaptive,dMat,verbose=FALSE) {
    bw <- bw.gwr3(formula,data=data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	gwr.model <- gwr.basic(formula,data=data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat)
	var.n <- length(gwr.model$lm$coefficients)
	coefs <- gwr.model$SDF@data[,1:var.n] # Key
	sds   <- gwr.model$SDF@data[,(var.n+6):(var.n*2+5)] # Key and corrected from CB & IG
	globcoef <- gwr.model$lm$coefficients # added by PH
	indep.vars <- colnames(gwr.model$lm$x)
	tvals <- c()
	for(i in 1:var.n)
	{
	  tvals <- cbind(tvals,(coefs[,i]-globcoef[i])/sds[,i]) # added by PH
	}
	#tvals <-as.data.frame(tvals)
	#names(tvals) <- paste(indep.vars, "MLR_t", sep="_")
	tvals
	}
gwrt.err <- function(data,formula,glw,approach, kernel, adaptive,dMat,verbose=FALSE) {
    bw <- bw.gwr3(formula,data=data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	gwr.model <- gwr.basic(formula,data=data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat)
	var.n <- length(names(gwr.model$lm$coefficients))
	coefs <- gwr.model$SDF@data[,1:var.n] # Key
	sds   <- gwr.model$SDF@data[,(var.n+6):(var.n*2+5)] # Key and corrected from CB & IG
	errmod <- errorsarlm(formula,data=data,listw=glw,method='spam')  # added by PH
	globcoef <- errmod$coefficients # added by PH
	indep.vars <- colnames(gwr.model$lm$x)
	tvals <- c()
	for(i in 1:var.n)
	{
	  tvals <- cbind(tvals,(coefs[,i]-globcoef[i])/sds[,i]) # added by PH
	}
	#tvals <-as.data.frame(tvals)
	#names(tvals) <- paste(indep.vars, "ERR_t", sep="_")
	tvals
	}
 #####
gwrt.lag <- function(data,formula,glw,approach, kernel, adaptive,dMat,verbose=FALSE) {
    bw <- bw.gwr3(formula,data=data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	gwr.model <- gwr.basic(formula,data=data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat)
	var.n <- length(names(gwr.model$lm$coefficients))
	coefs <- gwr.model$SDF@data[,1:var.n] # Key
	sds   <- gwr.model$SDF@data[,(var.n+6):(var.n*2+5)] # Key and corrected from CB & IG
	lagmod <- lagsarlm(formula,data=data,listw=glw,method='spam')  # added by PH
	globcoef <- lagmod$coefficients # added by PH
	indep.vars <- colnames(gwr.model$lm$x)
	tvals <- c()
	for(i in 1:var.n)
	{
	  tvals <- cbind(tvals,(coefs[,i]-globcoef[i])/sds[,i]) # added by PH
	}
	#tvals <-as.data.frame(tvals)
	#names(tvals) <- paste(indep.vars, "ERR_t", sep="_")
	tvals
	}

gwrt.sma <- function(data,formula, glw,approach, kernel, adaptive,dMat,verbose=FALSE) {
    bw <- bw.gwr3(formula,data=data,approach=approach,kernel=kernel, adaptive=adaptive,dMat=dMat,verbose=verbose)
	gwr.model <- gwr.basic(formula,data=data,bw=bw,kernel=kernel, adaptive=adaptive,dMat=dMat)
	var.n <- length(gwr.model$lm$coefficients)
	coefs <- gwr.model$SDF@data[,1:var.n] # Key
	sds   <- gwr.model$SDF@data[,(var.n+6):(var.n*2+5)] # Key and corrected from CB & IG
	olsmod <- lm(formula,data=data)
    smamod <- spautolm(formula,data=data,listw=glw,family='SMA')
    smamod$model <- olsmod$model
	globcoef <- summary(smamod)$Coef # added by PH
	indep.vars <- colnames(gwr.model$lm$x)
	tvals <- c()
	for(i in 1:var.n)
	{
	  tvals <- cbind(tvals,(coefs[,i]-globcoef[i,1])/sds[,i]) # added by PH
	}
	#tvals <-as.data.frame(tvals)
	#names(tvals) <- paste(indep.vars, "SMA_t", sep="_")
	tvals
	}
	
bw.gwr3<-function(formula, data, approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,verbose=FALSE,dMat, nlower = 10)
{
    ##Data points{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given regression data must be Spatial*DataFrame")
  }
  #cat("This selection has been optimised by golden selection.\n")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
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
  #################### Recommond to specify a distance matrix
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA
    if(approach=="cv"||approach=="CV")
       bw <- gold(gwr.cv,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat, dMat=dMat,verbose=verbose)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw<-gold(gwr.aic,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat,dMat=dMat,verbose=verbose)    
   # bw<-NA
#    if(approach=="cv"||approach=="CV")
#       bw <- optimize(bw.cv,lower=lower,upper=upper,maximum=FALSE,X=x,Y=y,kernel=kernel,
#       adaptive=adaptive, dp.locat=dp.locat, p=p, theta=theta, longlat=longlat,dMat=dMat,tol=.Machine$double.eps^0.25)
#    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
#       bw<-optimize(bw.aic,lower=lower,upper=upper,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
    bw

}