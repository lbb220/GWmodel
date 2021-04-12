## Don't need locations if dMat and dMat.rp are supplied
## Removed other inputs to gwr.mixed that weren't used in gwr.mixed
## Input diagnostic wasn't used in gwr.mixed, but is here:
##      if (diagnostic == F) this is now very fast
##      if (diagnostic == T) this is around 5 times faster than gwr.mixed
## Revised by Fiona H Evans from Centre for Digital Agriculture, Murdoch and Curtin Universities

gwr.mixed <- function(formula, data, regression.points, fixed.vars,intercept.fixed=FALSE, bw, diagnostic=T,
             kernel="bisquare", adaptive=FALSE, p=2, theta=0, longlat=F,dMat, dMat.rp)
{
   ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  if (diagnostic) 
     hatmatrix <- T
  else 
     hatmatrix <- F
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
  	rp.given <- FALSE
    regression.points <- data
    rp.locat<-coordinates(data)
  }
  else
  {
    rp.given <- TRUE
    if (is(regression.points, "Spatial"))
    {
       rp.locat<-coordinates(regression.points)
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 2)
       rp.locat<-regression.points
    else
      {
        warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
        rp.locat<-dp.locat
      }
  }
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
    #########Distance matrix is given or not
  dp.n <- nrow(dp.locat)
  rp.n <- nrow(rp.locat)
  if (missing(dMat))
  {
      DM.given<-F
      DM1.given<-F
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
      dMat.rp <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
       stop("Dimensions of dMat are not correct")
    if (missing(dMat.rp)) {
    dMat.rp <- dMat
    }
    else
    {
       dim.dMat.rp <- dim(dMat.rp)
    if (dim.dMat.rp[1]!=dp.n||dim.dMat.rp[2]!=rp.n)
        stop("Dimensions of dMat.rp are not correct")
    }
    DM1.given<-T 
  }
  ####################
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
    idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
      colnames(x)[idx1]<-"Intercept" 
  #colnames(x)[1]<-"Intercept"
  if (missing(fixed.vars))
  {
    warning("No independent variables in the formula is specified as fixed terms!")
    if(!intercept.fixed)
       stop("Please use basic GWR function to calibrate this model")
  }
  else
  {
     if(intercept.fixed)
        fixed.vars <- c("Intercept", fixed.vars)
  }
  idx.fixed <- match(fixed.vars, colnames(x))
  x1 <- x[, -idx.fixed]
  x2<- x[, idx.fixed]
  if (!is.null(x1)) x1 <- as.matrix(x1, nrow = dp.n)
  if (!is.null(x2)) x2 <- as.matrix(x2, nrow = dp.n)
  colnames(x1) <- colnames(x)[-idx.fixed]
  colnames(x2) <- colnames(x)[idx.fixed]
  #y <- as.matrix(y, nrow = dp.n)
  model <- gwr.mixed.2.fast(x1, x2, y, adaptive=adaptive, bw=bw,
                        kernel=kernel, dMat=dMat, dMat.rp=dMat.rp)                     
  res <- list()
   res$local <- model$local 
   res$global <- as.matrix(apply(model$global,2,mean,na.rm=T), 1, length(idx.fixed))
   colnames(res$local) <- colnames(x1)
   colnames(res$global) <- colnames(x2)
   mgwr.df <- data.frame(model$local, model$global)
   colnames(mgwr.df) <- c(paste(colnames(x1), "L", sep="_"), paste(colnames(x2), "F", sep="_"))
   rownames(rp.locat)<-rownames(mgwr.df)
  griddedObj <- F
     if (is(regression.points, "Spatial"))
     { 
         if (is(regression.points, "SpatialPolygonsDataFrame"))
         {
            polygons<-polygons(regression.points)
            #SpatialPolygons(regression.points)
   #         #rownames(mgwr.df) <- sapply(slot(polygons, "polygons"),
   #                            #  function(i) slot(i, "ID"))
            SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=mgwr.df, match.ID=F)
         }
         else
         {
            griddedObj <- gridded(regression.points)
            SDF <- SpatialPointsDataFrame(coords=rp.locat, data=mgwr.df, proj4string=CRS(p4s), match.ID=F)
            gridded(SDF) <- griddedObj 
         }
     }
     else
         SDF <- SpatialPointsDataFrame(coords=rp.locat, data=mgwr.df, proj4string=CRS(p4s), match.ID=F)
 # 
#   if (is(regression.points, "SpatialPolygonsDataFrame"))
#    {
#       polygons<-polygons(regression.points)
#       #SpatialPolygons(regression.points)
#       rownames(mgwr.df) <- sapply(slot(polygons, "polygons"),
#                          function(i) slot(i, "ID"))
#       SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=mgwr.df)
#    }
#    else
#       SDF <- SpatialPointsDataFrame(coords=rp.locat, data=mgwr.df, proj4string=CRS(p4s), match.ID=F)
   
   res$SDF <- SDF
   if (hatmatrix)
   {
      gwr.fitted <- function(x,b) apply(x*b,1,sum)
      edf <- gwr.mixed.trace.fast(x1, x2, y, adaptive=adaptive, bw=bw,
               kernel=kernel, dMat=dMat)
      model2 <-gwr.mixed.2.fast(x1, x2, y, adaptive=adaptive, bw=bw, 
                kernel=kernel, dMat=dMat, dMat.rp=dMat)
      #r.ss <- rss(y, cbind(x1,x2), cbind(model2$local, model2$global)) 
      r.ss <- sum((y - gwr.fitted(model2$global, x2) - gwr.fitted(model2$local,x1))^2)
      n1 <- length(y)
      sigma.aic <- r.ss / n1
      aic <- log(sigma.aic*2*pi) + 1 + 2*(edf + 1)/(n1 - edf - 2)
      aic <- n1*aic
      res$aic <- aic
      res$bic <- n1*log(sigma.aic) + n1*log(2*pi) + edf * log(n1)
      res$df.used <- edf
      res$r.ss <- r.ss
   }
   GW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,bw=bw, 
                       kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat,
                       DM.given=DM1.given,diagnostic=diagnostic)
   res$GW.arguments <- GW.arguments
   res$this.call <- this.call
   timings[["stop"]] <- Sys.time()
   res$timings <- timings
   class(res) <- "mgwr"
   res
}

gwr.mixed.2.fast <- function(x1, x2, y, adaptive=F, bw,
                            kernel="bisquare", dMat, dMat.rp)
{ 
  gwr_mixed_2(x1, x2, y, dMat, dMat.rp, bw, kernel, adaptive)
}


##Mixed GWR
gwr.mixed.2 <- function(x1, x2, y, loc, out.loc, adaptive=F, bw=sqrt(var(loc[,1])+var(loc[,2])),
               kernel="bisquare", p=2, theta=0, longlat=F,dMat)
{
  #gwr.fitted <- function(x,b) apply(x*b,1,sum)
  dp.n <- nrow(loc)
   
   ncols.2 <- dim(x2)[2]
   x3 <- NULL
   if(missing(out.loc))
     dMat1 <- dMat
   else
     dMat1 <- gw.dist(loc, p=p, theta = theta, longlat = longlat)
   for (i in 1:ncols.2)
   {
      m.temp <-gwr.q(x1, x2[,i], loc, adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat, dMat = dMat1)
      x3 <- cbind(x3,x2[,i]-gw.fitted(x1,m.temp))
   }
   colnames(x3) <- colnames(x2)
   m.temp <-gwr.q(x1, y, loc, adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat, dMat = dMat1)
   y2 <- y - gw.fitted(x1,m.temp)
   
   model2 <-gwr.q(x3, y2, loc, adaptive=TRUE, bw=1.0e6, kernel="boxcar",
                  p=p, theta=theta, longlat=longlat, dMat = dMat1)
   fit2 <- gw.fitted(x2,model2)
   model1 <-gwr.q(x1, y-fit2, loc, out.loc=out.loc,adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat,dMat=dMat)
   if(!missing(out.loc))
      model2 <-gwr.q(x3, y2, loc,out.loc=out.loc, adaptive=TRUE, bw=1.0e6, kernel="boxcar",
                  p=p, theta=theta, longlat=longlat,dMat=dMat)
   list(local=model1,global=model2)
  }

#
# Fix the column names for x3
#

gwr.mixed.trace <- function(x1, x2, y, loc, out.loc, adaptive=F, bw=sqrt(var(loc[,1])+var(loc[,2])),
               kernel="bisquare", p=2, theta=0, longlat=F,dMat)
  {gw.fitted <- gwr.fitted <- function(x,b) apply(x*b,1,sum)
   e.vec <- function(m,n) as.numeric(m == 1:n)
   dp.n <- nrow(loc)
   if (missing(dMat))
     DM.given <- F
   else
     DM.given <- T

   ncols.2 <- dim(x2)[2]
   #n.items <- length(y)
   if(missing(out.loc))
     dMat1 <- dMat
   else
     dMat1 <- gw.dist(loc, p=p, theta = theta, longlat = longlat)
   x3 <- NULL
   for (i in 1:ncols.2)
     {m.temp <-gwr.q(x1, x2[,i], loc, adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat,dMat=dMat1) 
      x3 <- cbind(x3,x2[,i]-gw.fitted(x1,m.temp))}
    
   colnames(x3) <- colnames(x2)
   hii <- NULL

   for (i in 1:dp.n)
     {
       m.temp <-gwr.q(x1, e.vec(i,dp.n), loc, adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat,dMat=dMat1)  
       y2 <- e.vec(i,dp.n) - gwr.fitted(x1,m.temp)

       model2 <-gwr.q(x3,y2, loc,  adaptive=TRUE, bw=1.0e6, kernel="boxcar",
                p=p, theta=theta, longlat=longlat,dMat=dMat1)
       fit2 <- gwr.fitted(x2,model2)
       if(DM.given)
       {
          model1 <-gwr.q(x1, e.vec(i,dp.n)-fit2, loc, out.loc=matrix(loc[i,], ncol=2), adaptive=adaptive, bw=bw,
                  kernel=kernel, dMat=matrix(dMat[,i], ncol=1))
          
          model2 <-gwr.q(x3,y2, loc, out.loc=matrix(loc[i,], ncol=2),  adaptive=TRUE, bw=1.0e6, kernel="boxcar",
                        dMat=matrix(dMat[,i], ncol=1))
          
       }
       else
       {
          model1 <-gwr.q(x1, e.vec(i,dp.n)-fit2, loc, out.loc=matrix(loc[i,], ncol=2), adaptive=adaptive, bw=bw,
                  kernel=kernel, p=p, theta=theta, longlat=longlat)
          model2 <-gwr.q(x3,y2, loc, out.loc=matrix(loc[i,], ncol=2),  adaptive=TRUE, bw=1.0e6, kernel="boxcar",
                        p=p, theta=theta, longlat=longlat)
       }  
       hii <- c(hii,gwr.fitted(matrix(x1[i,],nrow=1),model1)+gwr.fitted(matrix(x2[i,],nrow=1),model2)) }                   
   sum(hii)
  }
gwr.mixed.trace.fast <- function(x1, x2, y, adaptive=F, bw,
                            kernel="bisquare", dMat){
  gwr_mixed_trace(x1, x2, y, dMat, bw, kernel, adaptive)
}

print.mgwr <- function(x, ...)
{
  if(class(x) != "mgwr") stop("It's not a mgwr object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  
  cat("\n   *********************Model calibration information*********************\n")
  gwr.names <- colnames(x$local)
   global.names <- colnames(x$global)
   cat("   Mixed GWR model with local variables :", gwr.names, "\n")
   cat("   Global variables :", global.names, "\n")
	cat("   Kernel function:", x$GW.arguments$kernel, "\n")
	if(x$GW.arguments$adaptive)
	   cat("   Adaptive bandwidth: ", x$GW.arguments$bw, " (number of nearest neighbours)\n", sep="") 
  else
     cat("   Fixed bandwidth:", x$GW.arguments$bw, "\n")
	if(x$GW.arguments$rp.given) 
     cat("   Regression points: A seperate set of regression points is used.\n")
  else
     cat("   Regression points: the same locations as observations are used.\n")
	if (x$GW.arguments$DM.given) 
     cat("   Distance metric: A distance matrix is specified for this model calibration.\n")
  else
     {
     if (x$GW.arguments$longlat)
        cat("   Distance metric: Great Circle distance metric is used.\n")
     else if (x$GW.arguments$p==2)
        cat("   Distance metric: Euclidean distance metric is used.\n")
     else if (x$GW.arguments$p==1)
        cat("   Distance metric: Manhattan distance metric is used.\n") 
     else if (is.infinite(x$GW.arguments$p))
        cat("   Distance metric: Chebyshev distance metric is used.\n")
     else 
        cat("   Distance metric: A generalized Minkowski distance metric is used with p=",x$GW.arguments$p,".\n")
     if (x$GW.arguments$theta!=0&&x$GW.arguments$p!=2&&!x$GW.arguments$longlat)
        cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GW.arguments$theta, "in radian.\n")   
     }
   cat("\n   ****************Summary of mixed GWR coefficient estimates:******************\n")       
	 cat("   Estimated global variables :\n")
	 
	 gCM <- matrix(x$global,nrow=1)
	 rownames(gCM) <- "   Estimated global coefficients:" 
   colnames(gCM) <- global.names 
   printCoefmat(gCM)
   cat("   Estimated GWR variables :\n")
   	CM <- t(apply(x$local, 2, summary))[,c(1:3,5,6)]
   if(!is.matrix(CM))
   {
      CM <- as.matrix(t(CM), nrow=1)
      rownames(CM) <- gwr.names[1]
   }
    rnames<-rownames(CM)
    for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	 rownames(CM) <-rnames 
	 printCoefmat(CM)
   if (x$GW.arguments$diagnostic)
   {
    	cat("   ************************Diagnostic information*************************\n")
     cat("   Effective D.F.:  ",format(x$df.used,digits=4),"\n")
     cat("   Corrected AIC:  ",format(x$aic,digits=4),"\n")
     cat("             BIC:  ",format(x$bic,digits=4),"\n")
     cat("   Residual sum of squares:  ",format(x$r.ss,digits=4),"\n")
   }
  cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
   
}
   
gwr.q <- function(x, y, loc, out.loc=loc, adaptive=F, bw=sqrt(var(loc[,1])+var(loc[,2])),
                  kernel, p, theta, longlat,dMat, wt2=rep(1,nrow(loc)))
{
  if(missing(out.loc))
    rp.n <- nrow(loc)
  else
    rp.n <- nrow(out.loc)
  if (missing(dMat))
     dMat <- gw.dist(loc, out.loc, p, theta, longlat)
  var.n <- ncol(x)
  betas <- gwr_q(x,  y, dMat, bw, kernel, adaptive)
  colnames(betas) <- colnames(x)
  betas
}

gwr.q.fast <- function(x, y, adaptive=F, bw,
                  kernel, dMat)
{
  if (missing(dMat)) {
    stop("Error: distance matrix missing")
  }
  
  else {
    betas <- gwr_q(x,  y, dMat, bw, kernel, adaptive)
  }
    
  colnames(betas) <- colnames(x)
  betas
}