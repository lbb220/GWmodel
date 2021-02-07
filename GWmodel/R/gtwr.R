#Geographically and temporally weighted regression
###########################################
#Bo Huang , Bo Wu & Michael Barry (2010) Geographically and temporally weighted regression for modeling spatio-temporal 
#variation in house prices, International Journal of Geographical Information Science, 24:3, 383-401
# Bo Wu, Rongrong Li & Bo Huang (2014) A geographically and temporally
# weighted autoregressive model with application to housing prices, International Journal of Geographical Information Science, 28:5, 1186-1204

#Then an optimal spatial bandwidth is specified for each time period based on a goodnessof-
#  fit criterion such as cross-validation (CV) or Akaike information criterion (AIC). Using the
# optimal spatial bandwidth, the optimal temporal bandwidth is determined, again based on CV or
# AIC. Once both optimal spatial and temporal bandwidths are derived, they can be used to
# construct the spatiotemporal weight matrix W, which allows local parameters to be estimated
# using equation (2).

###########################################
#Calibrate the GTWR model
gtwr<- function(formula, data, regression.points, obs.tv, reg.tv, st.bw, kernel="bisquare",
                 adaptive=FALSE, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  polygons <- NULL
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    if(is(data, "SpatialPolygonsDataFrame"))
       polygons <- polygons(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
  dp.n <- nrow(dp.locat)
  ####Check the time stamps given for the data
  if(missing(obs.tv))
  {
    stop("Please provide the corresponding time stamps for the observations!")
  }
  else
  {
    if(is(obs.tv, "Date")||is(obs.tv, "POSIXlt")||is(obs.tv, "POSIXct")||is(obs.tv, "numeric")||is(obs.tv, "yearmon")||is(obs.tv, "yearqtr"))
    {
      if(length(obs.tv)!=dp.n)
        stop("The given time stamps must correspond strictly to the observation data")
    }
    else
    {
      stop("Please provide the time stamps in accepted format: numeric, Date, POSIXlt, POSIXct, yearmon or yearqtr")
    }
  }
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
    rp.given <- FALSE
    rp.locat<-dp.locat
    hatmatrix<-T
    reg.tv <- obs.tv
  }
  else 
  {
    rp.given <- TRUE
    hatmatrix<-F
    if (is(regression.points, "Spatial"))
    {
      rp.locat<-coordinates(regression.points)
      if (is(regression.points, "SpatialPolygonsDataFrame"))
         polygons<-polygons(regression.points)
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 2)
      rp.locat<-regression.points
    
    else
    {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      rp.given <- F
      rp.locat<-dp.locat
      reg.tv <- obs.tv
    }
    rp.n <- nrow(rp.locat)
    ###time stamps for regression locations
    if(missing(reg.tv))
      stop("Please provide the corresponding time stamps for the regression points!")
    else
    {
      if(is(reg.tv, "Date")||is(reg.tv, "POSIXlt")||is(reg.tv, "POSIXct")||is(reg.tv, "numeric")||is(reg.tv, "yearmon")||is(reg.tv, "yearqtr"))
      {
        if(length(reg.tv)!=rp.n)
          stop("The given time stamps must correspond strictly to the regression data")
      }
      else
      {
        stop("Please provide the time stamps in accepted format: numeric, Date, POSIXlt, POSIXct, yearmon or yearqtr")
      }
    }
  }
 #  if(class(obs.tv)==class(reg.tv))
 # {
 #   if(class(obs.tv)=="numeric" && t.units != "years")
#      t.units <- "auto"
#    if (class(obs.tv)=="yearmon")
#      t.units <- "months"
#    if (class(obs.tv)=="yearqtr")
#      t.units <- "quarters"
#  }
  ###
  ####################
  ######Extract the data frame
  ####Refer to the function lm
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n<-ncol(x)
  rp.n<-nrow(rp.locat)
  dp.n<-nrow(data)
  betas <-matrix(nrow=rp.n, ncol=var.n)
  betas.SE <-matrix(nrow=rp.n, ncol=var.n)
  betas.TV <-matrix(nrow=rp.n, ncol=var.n)
  ##S: hatmatrix
  S<-matrix(nrow=dp.n,ncol=dp.n)
  #C.M<-matrix(nrow=dp.n,ncol=dp.n)
  idx1 <- match("(Intercept)", colnames(x))
  if(!is.na(idx1))
    colnames(x)[idx1]<-"Intercept" 
  colnames(betas) <- colnames(x)
  #colnames(betas)[1]<-"Intercept"
  
  ##################################################
  #####Linear regression
  lm.res <- lm(formula,data=data)
  lm.res$x <- x
  lm.res$y <- y
  gTSS <- c(cov.wt(matrix(y, ncol=1), wt=rep(as.numeric(1), dp.n), method="ML")$cov*dp.n)
  #####GTWR
  ######### Spatial Distance matrix is given or not
  if (missing(st.dMat))
  {
    DM.given<-F
    if(dp.n+rp.n <10000)
    {
      if(rp.given)
        st.dMat <- st.dist(dp.locat, rp.locat, obs.tv, reg.tv, p=p, theta=theta, longlat=longlat,lamda=lamda,t.units = t.units,ksi=ksi)
      else
        st.dMat <- st.dist(dp.locat, obs.tv=obs.tv, p=p, theta=theta, longlat=longlat,lamda=lamda,t.units = t.units,ksi=ksi)
      DM.given <- T
    }
  }
  else
  {
    DM.given<-T
    dim.stdMat<-dim(st.dMat)
    if (dim.stdMat[1]!=dp.n||dim.stdMat[2]!=rp.n)
      stop("Dimensions of spatio-temporal distance matrix sdMat are not correct")
  }
  #############Calibration the model
  for (i in 1:rp.n)
  {
    if(DM.given)
      st.disti <- st.dMat[,i]
    else
      st.disti <- st.dist(dp.locat, rp.locat, obs.tv, reg.tv, focus=i,p=p, theta=theta, longlat=F,lamda=longlat,t.units = t.units,ksi=ksi)
    W.i<-gw.weight(st.disti,st.bw,kernel,adaptive)
    gw.resi<-gw_reg(x,y,W.i,hatmatrix,i)
    betas[i,]<-gw.resi[[1]] ######See function by IG
    if(hatmatrix)
    {
      S[i,]<-gw.resi[[2]]
      Ci<-gw.resi[[3]]
      betas.SE[i,]<-diag(Ci%*%t(Ci))
    }
  }
  ########################Diagnostic information
  
  GTW.diagnostic<-NA
  if (hatmatrix)
  {
    tr.S<-sum(diag(S))
    tr.StS<-sum(S^2)
    Q<-t(diag(dp.n)-S)%*%(diag(dp.n)-S)
    RSS.gw<-t(y)%*%Q%*%y
    yhat<-S%*%y
    residual<-y-yhat
    #####Calculate the standard errors of the parameter estimates
    #pseudo-t values
    sigma.hat1<-RSS.gw/(dp.n-2*tr.S+tr.StS)
    Stud_residual<-residual
    q.diag<-diag(Q)
    for(i in 1:dp.n)
    {
      Stud_residual[i]<-residual[i]/sqrt(sigma.hat1*q.diag[i])
      betas.SE[i,]<-sqrt(sigma.hat1*betas.SE[i,])
      betas.TV[i,]<-betas[i,]/betas.SE[i,] 
    }
    sigma.hat2 <- RSS.gw/dp.n
    AIC<-dp.n*log(sigma.hat2) + dp.n*log(2*pi) +dp.n+tr.S
    AICc<-dp.n*log(sigma.hat2) + dp.n*log(2*pi) + dp.n *((dp.n + tr.S) / (dp.n - 2 - tr.S))
    edf<- dp.n - 2*tr.S + tr.StS
    enp<-2*tr.S - tr.StS
    yss.g <- sum((y - mean(y))^2)
    gw.R2<-1-RSS.gw/yss.g; ##R Square valeu
    gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared value
    GTW.diagnostic<-list(RSS.gw=RSS.gw,AIC=AIC,AICc=AICc,enp=enp, edf=edf,gw.R2=gw.R2,gwR2.adj=gwR2.adj)
  }
  
  ####encapsulate the GWR results
  GTW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,st.bw=st.bw,lamda=lamda, ksi=ksi,kernel=kernel,
                      adaptive=adaptive,p=p, theta=theta, longlat=longlat,DM.given=DM.given,units=units)
  
  #observed y: y
  #fitted y : yhat 
  #residual : y-yhat
  #Studentised residual: Stud_residual=residual.i/(sigma.hat*sqrt(q.ii))
  if (hatmatrix)                                         
  {
    gtwres.df<-data.frame(betas,y,yhat,residual,obs.tv, Stud_residual,betas.SE,betas.TV)
    colnames(gtwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual","time_stamp","Stud_residual")),
                             paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    
  }
  else
  {
    
    gtwres.df<-data.frame(betas, reg.tv)
    colnames(gtwres.df)<- c(colnames(betas),"time_stamp") 
  }
  rownames(rp.locat)<-rownames(gtwres.df)
  

  if (!is.null(polygons))
  {
    rownames(gtwres.df) <- sapply(slot(polygons, "polygons"),
                                  function(i) slot(i, "ID"))
    SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gtwres.df,match.ID=F)
  }
  else
  {
    SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gtwres.df, proj4string=CRS(p4s), match.ID=F)
  }
  timings[["stop"]] <- Sys.time()
  ##############
  res<-list(GTW.arguments=GTW.arguments,GTW.diagnostic=GTW.diagnostic,lm=lm.res,SDF=SDF,
            timings=timings,this.call=this.call)
  class(res) <-"gtwrm"
  invisible(res) 
}

############################Layout function for outputing the GWR results
##Author: BL	
print.gtwrm<-function(x, ...)
{
  if(class(x) != "gtwrm") stop("It's not a gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GTW.arguments$formula)
  var.n<-length(x$lm$coefficients)
  cat("\n   Dependent (y) variable: ",vars[1])
  cat("\n   Independent variables: ",vars[-1])
  dp.n<-length(x$lm$residuals)
  cat("\n   Number of data points:",dp.n)
  ################################################################ Print Linear
  cat("\n   ***********************************************************************\n")
  cat("   *                    Results of Global Regression                     *\n")
  cat("   ***********************************************************************\n")
  print(summary.lm(x$lm))
  cat("   ***Extra Diagnostic information\n")
  lm_RSS<-sum(x$lm$residuals^2)
  lm_Rank<-x$lm$rank     
  cat("   Residual sum of squares:", lm_RSS)
  #lm_sigma<-sqrt(lm_RSS/(dp.n-lm_Rank-2))
  lm_sigma<-sqrt(lm_RSS/(dp.n-2))
  cat("\n   Sigma(hat):", lm_sigma)
  lm_AIC<-dp.n*log(lm_RSS/dp.n)+dp.n*log(2*pi)+dp.n+2*(var.n + 1)
  #AIC = dev + 2.0 * (double)(MGlobal + 1.0);
  cat("\n   AIC: ", lm_AIC)
  ##AICc = 	dev + 2.0 * (double)N * ( (double)MGlobal + 1.0) / ((double)N - (double)MGlobal - 2.0);
  lm_AICc= dp.n*log(lm_RSS/dp.n)+dp.n*log(2*pi)+dp.n+2*dp.n*(var.n+1)/(dp.n-var.n-2)
  cat("\n   AICc: ", lm_AICc)
  #lm_rdf <- x$dfsidual
  
  #########################################################################
  cat("\n   ***********************************************************************\n")
    cat("   *    Results of Geographically and Temporally Weighted Regression     *\n")
  cat("   ***********************************************************************\n")
  cat("\n   *********************Model calibration information*********************\n")
  cat("   Kernel function for geographically and temporally weighting:", x$GTW.arguments$kernel, "\n")
  if(x$GTW.arguments$adaptive)
    cat("   Adaptive bandwidth for geographically and temporally  weighting: ", x$GTW.arguments$st.bw, " (number of nearest neighbours)\n", sep="")
  else
    cat("   Fixed bandwidth for geographically and temporally weighting: ", x$GTW.arguments$st.bw, "\n")
  if(x$GTW.arguments$rp.given) 
    cat("   Regression points: A seperate set of regression points is used.\n")
  else
    cat("   Regression points: the same locations as observations are used.\n")
  if (x$GTW.arguments$DM.given)
    cat("   Distance metric for geographically and temporally  weighting: A distance matrix is specified for this model calibration.\n")
  else
  {
    if (x$GTW.arguments$longlat)
      cat("   Distance metric for geographically weighting: Great Circle distance metric is used.\n")
    else if (x$GTW.arguments$p==2)
      cat("   Distance metric for geographically weighting: Euclidean distance metric is used.\n")
    else if (x$GTW.arguments$p==1)
      cat("   Distance metric for geographically weighting: Manhattan distance metric is used.\n") 
    else if (is.infinite(x$GTW.arguments$p))
      cat("   Distance metric for geographically weighting: Chebyshev distance metric is used.\n")
    else 
      cat("   Distance metric for geographically weighting: A generalized Minkowski distance metric is used with p=",x$GTW.arguments$p,".\n")
    if (x$GTW.arguments$theta!=0&&x$GTW.arguments$p!=2&&!x$GTW.arguments$longlat)
      cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GTW.arguments$theta, "in radian.\n")
    cat("   The temporal distance is calculated in ", x$GTW.arguments$units, ".\n")
    cat("   The adjustment parameter for calculating spatio-temporal distances lamda is: ", x$GTW.arguments$lamda, ".\n")
    cat("   The adjustment parameter for calculating spatio-temporal distances ksi is: ", x$GTW.arguments$ksi, "in radian.\n")    
  } 
  
  cat("\n   ****************Summary of GTWR coefficient estimates:*****************\n")       
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
  if (x$GTW.arguments$hatmatrix) 
  {	
    cat("   ************************Diagnostic information*************************\n")
    cat("   Number of data points:", dp.n, "\n")
    cat("   Effective number of parameters (2trace(S) - trace(S'S)):", x$GTW.diagnostic$enp, "\n")
    cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GTW.diagnostic$edf, "\n")
    cat("   AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33):",
        x$GTW.diagnostic$AICc, "\n")
    cat("   AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22):", x$GTW.diagnostic$AIC, "\n")
    cat("   Residual sum of squares:", x$GTW.diagnostic$RSS.gw, "\n")
    cat("   R-square value: ",x$GTW.diagnostic$gw.R2,"\n")
    cat("   Adjusted R-square value: ",x$GTW.diagnostic$gwR2.adj,"\n")	
  }
  cat("\n   ***********************************************************************\n")
  cat("   Program stops at:", as.character(x$timings$stop), "\n")
  invisible(x)
}

#Calculate the time distance vector 
ti.distv <- function(focal.t, obs.tv, units="auto")
{
  n <- length(obs.tv)
  dist.tv <- c()
  for (t in obs.tv) 
  {
    if(focal.t >=t)
       dist.tv <- c(dist.tv, ti.dist(t,focal.t,units = units))
    else
       dist.tv <- c(dist.tv, Inf)
  }
  dist.tv
}
#Calculate the time distance matrix
ti.distm <- function(obs.tv,reg.tv, units="auto")
{
  n <- length(obs.tv)
  if(missing(reg.tv))
  {
    m.sys <- T
    m <- n
    reg.tv <- obs.tv
  }
  else
  {
    m.sys <- F
    m <- length(reg.tv)
  }
  dist.tm <- matrix(numeric(m*n),nrow=n)
  #if(m.sys)
#  {
#    for(i in 1:n)
#      for(j in 1:i)
#      {
#        dist.tm[i,j] <- ti.dist(obs.tv[i],obs.tv[j],units = units)
#        dist.tm[j,i] <- dist.tm[i,j]
#      }
#  }
#  else
#  {
#    for (i in 1:m) 
#      for (j in 1:n) 
#      {
#        dist.tm[i,j] <- ti.dist(obs.tv[i],obs.tv[j],units = units)
#      }
#  }
  for(i in 1:m)
    dist.tm[,i] <- ti.distv(reg.tv[i], obs.tv, units)
  dist.tm
}
#calculate the time distance
#units can be "auto", "secs", "mins", "hours","days", "weeks","months","years"
ti.dist <- function(t1,t2,units="auto")
{
  tcl <- class(t1)
  switch(tcl,
         Date = as.numeric(difftime(t1,t2,units = units)),
         POSIXlt = as.numeric(difftime(t1,t2,units = units)),
         POSIXct = as.numeric(difftime(t1,t2,units = units)),
         numeric = t2-t1,
         integer = t2-t1,
         yearmon = (t2-t1)*12,
         yearqtr = (t2-t1)*4
         )
}
#Get the forward time stamps before a specific time stamp
get.before.ti <- function(ti, ts, time.lag,units)
{
  ts.ti <- c()
  idx <- which(ts==ti)
  idx <- idx-1
  while(idx>0 && ti.dist(ti, ts[idx], units)<time.lag)
  {
    ts.ti <- c(ts.ti,ts[idx])
    idx <- idx-1
  }
  ts.ti
}
#Get the time stamps for the ST data frame, return the sorted time stamps and index
get.ts <- function(tv)
{
  t.fac <- as.factor(tv)
  ts <- as(sort(levels(t.fac)), class(tv))
  n <- length(ts)
  idxs <- c()
  for(i in tv)
  {
    idxs <- c(idxs, which(ts==i))
  }
  res <- list(ts=ts, index=idxs) 
  res
}
#Get the unique locations and index
get.uloat <- function(coords)
{
  TF.dup <- duplicated(coords)
  ucoords <- coords[!TF.dup,]
  n <- nrow(coords)
  idx <- c()
  for(i in 1:n)
  {
    idx <-c(idx, which(ucoords[,1]==coords[i,1]&ucoords[,2]==coords[i,2]))
  }
  res <- list(ucoords, idx)
  res
}
##Calculate the spatial distance matrix with duplicated locations removed
sdist.mat <- function(dp.locat, rp.locat, p=2, theta=0, longlat=F)
{
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
   {
       rp.given<-T
   }
   else
   {
       rp.given<-F 
       rp.locat<- dp.locat
   } 
   if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
   else
      rp.locat <- matrix(rp.locat, ncol=2)
   ###Remove the duplicated locations
   ucoord.dp <- get.uloat(dp.locat) 
   coord.dp.idx <- ucoord.dp[[2]]
   ucoord.dp <- ucoord.dp[[1]]
   ucoord.rp <- get.uloat(rp.locat) 
   coord.rp.idx <- ucoord.dp[[2]]
   ucoord.rp <- ucoord.rp[[1]]
   if(rp.given)
   {
     dists <- gw.dist(ucoord.dp, ucoord.rp, p=p, theta=theta, longlat=longlat) 
   }
   else
   {
     dists <- gw.dist(ucoord.dp, p=p, theta=theta, longlat=longlat) 
   }
   dists
}
tdist.mat <- function(obs.tv, reg.tv, t.units = "auto")
{
   uts.obv <- get.ts(obs.tv)
   uts.obv.idx <- uts.obv[[2]]
   uts.obv <- uts.obv[[1]]
   if(missing(reg.tv))
      dists <- ti.distm(uts.obv, units=t.units)
   else
   {
     uts.reg <- get.ts(reg.tv)
     uts.reg.idx <- uts.reg[[2]]
     uts.reg <- uts.reg[[1]]
     dists <- ti.distm(uts.obv,uts.reg, units=t.units)
   }
   dists
}
#Calculate the ST distance matrix
st.dist <- function(dp.locat, rp.locat, obs.tv, reg.tv,focus=0, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, s.dMat,t.dMat)
{
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
   {
       rp.given<-T
   }
   else
   {
       rp.given<-F 
       rp.locat<- dp.locat
   } 
   if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
   else
      rp.locat <- matrix(rp.locat, ncol=2)
   if (focus<0||focus>length(rp.locat[,1]))
      stop("No regression point is fixed")

   n.rp<-length(rp.locat[,1])
   n.dp<-length(dp.locat[,1])
   
   if (focus>0)
       dists<-numeric(n.dp) 
   else
       dists<-matrix(numeric(n.rp*n.dp),nrow=n.dp)
   ###Remove the duplicated locations
   ucoord.dp <- get.uloat(dp.locat) 
   coord.dp.idx <- ucoord.dp[[2]]
   ucoord.dp <- ucoord.dp[[1]]
   ucoord.rp <- get.uloat(rp.locat) 
   coord.rp.idx <- ucoord.rp[[2]]
   ucoord.rp <- ucoord.rp[[1]]
   ######Calculate the spatial distance matrix
   if(missing(s.dMat))
   {
      if(rp.given)
        s.dMat <- gw.dist(ucoord.dp, ucoord.rp, p=p, theta=theta, longlat=longlat)
      else
        s.dMat <- gw.dist(ucoord.dp, p=p, theta=theta, longlat=longlat) 
   }
   else
   {
      if(!(dim(s.dMat)[1]==nrow(ucoord.dp)&&dim(s.dMat)[2]==nrow(ucoord.rp)))
        stop("s.dMat is of dimnensions with duplicated locations removed")
   }
   
   ####Calculate the temporal distance matrix
   uts.obv <- get.ts(obs.tv)
   uts.obv.idx <- uts.obv[[2]]
   uts.obv <- uts.obv[[1]]
   if(rp.given)
   {
     uts.reg <- get.ts(reg.tv)
     uts.reg.idx <- uts.reg[[2]]
     uts.reg <- uts.reg[[1]]
   }
   if(missing(t.dMat))
   {
      if(rp.given)
        t.dMat <- ti.distm(uts.obv,uts.reg, units=t.units)
      else
        t.dMat <- ti.distm(uts.obv, units=t.units)
   }
   ####Calculate the distance matrix
   if(focus>0)
   {
     if(rp.given)
     {
       for(i in 1:n.dp)
       {
       if(is.infinite(t.dMat[uts.obv.idx[i], uts.reg.idx[focus]]))
            {
               dists[i] <- Inf
            }
       else
       {
         dists[i] <- lamda*s.dMat[coord.dp.idx[i],coord.rp.idx[focus]]+(1-lamda)*t.dMat[uts.obv.idx[i], uts.reg.idx[focus]]+2*sqrt(lamda*(1-lamda)*s.dMat[coord.dp.idx[i],coord.rp.idx[focus]]*t.dMat[uts.obv.idx[i], uts.reg.idx[focus]])*cos(ksi)
         }
       }
     }
     else
     {
       for(i in 1:n.dp)
       {
          if(is.infinite(t.dMat[uts.obv.idx[i], uts.obv.idx[focus]]))
            {
               dists[i] <- Inf
            }
            else
            {
         dists[i] <- lamda*s.dMat[coord.dp.idx[i],coord.dp.idx[focus]]+(1-lamda)*t.dMat[uts.obv.idx[i], uts.obv.idx[focus]]+2*sqrt(lamda*(1-lamda)*s.dMat[coord.dp.idx[i],coord.dp.idx[focus]]*t.dMat[uts.obv.idx[i], uts.obv.idx[focus]])*cos(ksi)
         }
       }
     }
   }
   else
   {
      if(rp.given)
      {
        for(j in 1:n.rp)
          for(i in 1:n.dp)
          {
            if(is.infinite(t.dMat[uts.obv.idx[i], uts.reg.idx[j]]))
            {
               dists[i,j] <- Inf
            }
            else
            {
            dists[i,j] <- lamda*s.dMat[coord.dp.idx[i],coord.rp.idx[j]]+(1-lamda)*t.dMat[uts.obv.idx[i], uts.reg.idx[j]]+2*sqrt(lamda*(1-lamda)*s.dMat[coord.dp.idx[i],coord.rp.idx[j]]*t.dMat[uts.obv.idx[i], uts.reg.idx[j]])*cos(ksi)
            }
          }
      }
     else
     {
       for(j in 1:n.dp)
         for(i in 1:j)
         {
           if(is.infinite(t.dMat[uts.obv.idx[i], uts.obv.idx[j]]))
            {
               dists[i,j] <- Inf
            }
            else
            {
           dists[i,j] <- lamda*s.dMat[coord.dp.idx[i],coord.dp.idx[j]]+(1-lamda)*t.dMat[uts.obv.idx[i], uts.obv.idx[j]]+2*sqrt(lamda*(1-lamda)*s.dMat[coord.dp.idx[i],coord.dp.idx[j]]*t.dMat[uts.obv.idx[i], uts.obv.idx[j]])*cos(ksi)
            }
           dists[j,i] <- dists[i,j]
         }
     }
   }
   dists
}

