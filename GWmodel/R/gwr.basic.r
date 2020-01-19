#########################
###Basic functions for GWR
# Author: BL, edited IG


##A function for calibrate a basic GWR model
#formula: Regression model formula of a formula object
#data: a dataframe, or a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp
#regression.points:an object containing the coordinates of regression points, often an object from package sp; if missing,the given data points are used as regression points
#bw: bandwidth used in the weighting function, possibly calculated by bw.sel()
####related with the function gw.weight
################################
##kernel: text vector of function chosen as follows
#	"gaussian": wgt = exp(-.5*(vdist/bw)^2)
# "bisquare": wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise
#	"tricube": wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise
# "boxcar": wgt=1 if dist < bw, wgt=0 otherwise
#adaptive:if TRUE calulate the adaptive kernel, and bw correspond to the number of nearest neighbours, default is FALSE.
#####################related with the function gw.dist
#p: the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
# theta: an angle in radian to roate the coordinate system, default is 0
# longlat: if TRUE, great circle will be caculated
# dMat: if given, it will be used for the GWR calibration

########################################################To be worked
#The outputs for each focus point should include:
#
#
#observed y
#fitted y
#residual
#standardised residual
#
#
#parameter estimates
#standard errors of the parameter estimates
#pseudo-t values
#
#Belsey-Kuh-Welsh condition number
#Variance Inflation Factors
#Variance decomposition proportions
gwr.basic <- function(formula, data, regression.points, bw, kernel="bisquare", adaptive=FALSE, p=2, theta=0, longlat=F, dMat, F123.test=F, cv=F, W.vect=NULL)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
  	rp.given <- FALSE
    regression.points <- data
    rp.locat <- coordinates(data)
    hatmatrix <- T
  }
  else
  {
    rp.given <- TRUE
    hatmatrix<-F
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
  #Regression data is gridded or not
  griddedObj <- F
  if(is(regression.points, "Spatial")) {
    if (is(regression.points, "SpatialPolygonsDataFrame")) polygons<-polygons(regression.points)
    else griddedObj <- gridded(regression.points)
  }
  ##Data points{
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  } else {
    stop("Given regression data must be Spatial*DataFrame")
  }

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
  if (hatmatrix) {
    betas.SE <-matrix(nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(nrow=rp.n, ncol=var.n)
  }
  idx1 <- match("(Intercept)", colnames(x))
  if(!is.na(idx1))
    colnames(x)[idx1]<-"Intercept"
  colnames(betas) <- colnames(x)
  # colnames(betas)[1]<-"Intercept"

  ##################################################
  ##### Linear regression
  lms <- lm(formula,data=data)
  # lms <-fastLm(formula,data=data)
  lms$x <- x
  lms$y <- y
  gTSS <- c(cov.wt(matrix(y, ncol=1), wt=rep(as.numeric(1), dp.n), method="ML")$cov*dp.n)

  ####################################################GWR
	  #########Distance matrix is given or not

  DM.given <- F
  if (missing(dMat))
  {
    DM.given<-F
    DM1.given<-F
    if(dp.n + rp.n <= 5000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    } else {
      dMat <- matrix(0, 1, 1)
    }
  }
  else
  {
    DM.given<-T
    DM1.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
       stop("Dimensions of dMat are not correct")
  }
  #############Calibration the model
  # W <- matrix(nrow = dp.n, ncol = rp.n)
  s_hat <- c(0.0, 0.0)
  q.diag <- matrix(0, 1, dp.n)
  for (i in 1:rp.n) {
    if (DM.given) dist.vi<-dMat[,i] else {
      if (rp.given) dist.vi<- gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat)
      else dist.vi<- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }
	W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    if (!is.null(W.vect)) W.i <- W.i * W.vect
    gwsi <- gw_reg(x, y, W.i, hatmatrix, i)
    betas[i,] <- gwsi[[1]] ######See function by IG
    if(hatmatrix)
    {
      si <- gwsi[[2]]
      Ci <- gwsi[[3]]
      betas.SE[i,] <- rowSums(Ci * Ci)
      s_hat[1] = s_hat[1] + si[i]
      s_hat[2] = s_hat[2] + sum(si %*% t(si))
      onei <- numeric(rp.n)
      onei[i] = 1
      p_i = onei - si
      q.diag = q.diag + p_i * p_i
    }
  }
  ########################Diagnostic information
  GW.diagnostic <- NA
  Ftests <- list()
  if (hatmatrix)
  {
    diags <- gwr_diag1(y, x, betas, as.vector(s_hat))
    tr.S <- s_hat[1]
    tr.StS <- s_hat[2]
    #print(paste("tr.S", tr.S))
    #print(paste("tr.StS", tr.StS))
    RSS.gw <- diags[5]
    yhat <- gw.fitted(x, betas)
    residual <- y - yhat
    CV <- numeric(dp.n)
    local.R2 <- numeric(dp.n)
    if (cv)
       CV<-gwr.cv.contrib(bw, x, y, kernel,adaptive, dp.locat, p, theta, longlat,dMat)

    #Studentised residual: Stud_residual=residual.i/(sigma.hat*sqrt(q.ii))
    #sigma.hat1=RSS/(n-2trace(S) + trace(S'S):Effective degrees of freedom
    #####Calculate the standard errors of the parameter estimates
    #pseudo-t values
    sigma.hat1 <- RSS.gw / (dp.n - 2 * tr.S + tr.StS)
    Stud_residual <- residual / t(sqrt(sigma.hat1 * q.diag))
    betas.SE <- sqrt(sigma.hat1 * betas.SE)
    betas.TV <- betas / betas.SE
    dybar2 <- (y - mean(y))^2
    dyhat2 <- (y - yhat)^2
    if (DM.given) {
      W <- gw.weight(dMat, bw, kernel, adaptive)
      TSSw <- W %*% dybar2
      RSSw <- W %*% dyhat2
      local.R2 <- (TSSw - RSSw) / TSSw
    } else {
      dybar2 <- t(dybar2)
      dyhat2 <- t(dyhat2)
      for(i in 1:dp.n)
      {
        if (DM.given) dist.vi <- dMat[,i] else {
          if (rp.given) dist.vi <- gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat)
          else dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        if (!is.null(W.vect)) W.i<-W.i*W.vect
        TSSw <- dybar2 %*% W.i
        RSSw <- dyhat2 %*% W.i
        local.R2[i] <- (TSSw - RSSw) / TSSw
      }
    }
    AIC <- diags[1]
    AICc <- diags[2]
    edf <- diags[3]
    enp <- diags[4]
    gw.R2 <- diags[6]
    gwR2.adj <- diags[7]
    GW.diagnostic <- list(RSS.gw = RSS.gw, AIC = AIC, AICc = AICc, enp = enp, edf = edf, gw.R2 = gw.R2, gwR2.adj = gwR2.adj)
    ######Parameters returned for F tests
    Ftests<-list()
    if(F123.test)
    {
      F.test.parameters <- list(
        dp.n = dp.n,
        var.n = var.n,
        dMat = dMat,
        dp.locat = dp.locat,
        x = x,
        bw = bw,
        adaptive = adaptive,
        kernel = kernel,
        betas = betas,
        RSS.lm = sum(lms$residuals^2),
        DF.lm = lms$df.residual,
        RSS.gw = RSS.gw,
        tr.S = tr.S,
        tr.StS = tr.StS,
        q.diag = q.diag,
        W.vect = W.vect,
        p = p,
        theta = theta,
        longlat = longlat
      )
      Ftests <- F1234.test(F.test.parameters)
    }
  }

  ####encapsulate the GWR results
  GW.arguments <- list(formula = formula, rp.given = rp.given, hatmatrix = hatmatrix, bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, theta = theta, longlat = longlat,
                       DM.given = DM1.given, F123.test = F123.test)

  #observed y: y
  #fitted y : yhat
  #residual : y-yhat
  #Studentised residual: Stud_residual=residual.i/(sigma.hat*sqrt(q.ii))
  if (hatmatrix)
  {
    if (is.null(W.vect))
    {
       gwres.df<-data.frame(betas,y,yhat,residual,CV,Stud_residual,betas.SE,betas.TV,local.R2)
       colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual","CV_Score","Stud_residual")),
                           paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"), "Local_R2")
    }
    else
    {
       gwres.df<-data.frame(betas,y,yhat,residual,CV,Stud_residual,betas.SE,betas.TV, W.vect,local.R2)
       colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual","CV_Score","Stud_residual")),
                          paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"), "E_weigts","Local_R2")
    }

  }
  else
  {
    if (is.null(W.vect))
        gwres.df<-data.frame(betas)
    else
    {
        gwres.df<-data.frame(betas, W.vect)
        colnames(gwres.df)<- c(colnames(betas), "E_weigts")
    }
  }
  rownames(rp.locat)<-rownames(gwres.df)

  if (is(regression.points, "SpatialPolygonsDataFrame"))
  {
     polygons<-polygons(regression.points)
     #SpatialPolygons(regression.points)
     rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                          function(i) slot(i, "ID"))
     SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df,match.ID=F)
  }
  else
  {
     SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
     if(griddedObj)
        gridded(SDF) <- T
  }
  timings[["stop"]] <- Sys.time()
 ##############
  res<-list(GW.arguments=GW.arguments,GW.diagnostic=GW.diagnostic,lm=lms,SDF=SDF,
            timings=timings,this.call=this.call,Ftests=Ftests)
  class(res) <-"gwrm"
  invisible(res)
}

reg.combine <- function(before, item) {
  result.names <- names(before)
  i <- item$index
  before$betas[i,] <- item$betas[i,]
  if (all(c("betas.SE", "s_hat", "q.diag") %in% result.names)) {
    before$betas.SE[i,] <- item$betas.SE[i,]
    before$s_hat = before$s_hat + item$s_hat
    before$q.diag = before$q.diag + item$q.diag
  }
  before
}


############################Layout function for outputing the GWR results
##Author: BL
print.gwrm<-function(x, ...)
{
  if(class(x) != "gwrm") stop("It's not a gwm object")
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
	  cat("   *          Results of Geographically Weighted Regression              *\n")
	cat("   ***********************************************************************\n")
	cat("\n   *********************Model calibration information*********************\n")
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
	if (x$GW.arguments$hatmatrix)
  {
    cat("   ************************Diagnostic information*************************\n")
		cat("   Number of data points:", dp.n, "\n")
		cat("   Effective number of parameters (2trace(S) - trace(S'S)):", x$GW.diagnostic$enp, "\n")
		cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GW.diagnostic$edf, "\n")
		cat("   AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33):",
                    x$GW.diagnostic$AICc, "\n")
		cat("   AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22):", x$GW.diagnostic$AIC, "\n")
		cat("   Residual sum of squares:", x$GW.diagnostic$RSS.gw, "\n")
    cat("   R-square value: ",x$GW.diagnostic$gw.R2,"\n")
		cat("   Adjusted R-square value: ",x$GW.diagnostic$gwR2.adj,"\n")
  }
  if(x$GW.arguments$F123.test && x$GW.arguments$hatmatrix)
		{
			cat("   ******************F test results of GWR calibration********************\n")
			cat("   ---F1 test (Leung et al. 2000)\n")
      rownames(x$Ftests$F1.test)<-"   "
			printCoefmat(x$Ftests$F1.test,signif.legend=F)
			cat("   ---F2 test (Leung et al. 2000)\n")
			rownames(x$Ftests$F2.test)<-"   "
			printCoefmat(x$Ftests$F2.test,signif.legend=F)
			cat("   ---F3 test (Leung et al. 2000)\n")
			rownames(x$Ftests$F3.test)<-rnames
			printCoefmat(x$Ftests$F3.test,signif.legend=F)
			cat("   ---F4 test (GWR book p92)\n")
			rownames(x$Ftests$F4.test)<-"   "
			printCoefmat(x$Ftests$F4.test,signif.legend=F)
			cat("\n   ---Significance stars")
			cat("\n   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ")
		}
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}
#######################################Re-organised the output of lm function
####Edited print the lm results
##Editor: BL
print.summary.lm <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
	      signif.stars= getOption("show.signif.stars"),	...)
{
    cat("\n   Call:\n", # S has ' ' instead of '\n'
	paste("   ",paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep=""))
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if(!is.null(x$w) && diff(range(x$w))) "Weighted ",
        "   Residuals:\n", sep="")
    if (rdf > 5L) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2L)
	    structure(apply(t(resid), 1L, quantile),
		      dimnames = list(nam, dimnames(resid)[[2L]]))
	else  {
            zz <- zapsmall(quantile(resid), digits + 1)
            structure(zz, names = nam)
        }
	print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("   ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
    }
    if (length(x$aliased) == 0L) {
        cat("\n   No Coefficients\n")
    } else {
        if (nsingular <- df[3L] - df[1L])
            cat("\n   Coefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\n   Coefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        rnames<-rownames(coefs)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	   rownames(coefs)<-rnames
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars, na.print="NA", signif.legend=F)
		cat("\n   ---Significance stars")
		cat("\n   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ")
    }
    ##
    cat("\n   Residual standard error:",
	format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    if (!is.null(x$fstatistic)) {
	cat("   Multiple R-squared:", formatC(x$r.squared, digits=digits))
	cat("\n   Adjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
	    "\n   F-statistic:", formatC(x$fstatistic[1L], digits=digits),
	    "on", x$fstatistic[2L], "and",
	    x$fstatistic[3L], "DF,  p-value:",
	    format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE), digits=digits),
	    "\n")
    }
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1L) {
	    cat("\n   Correlation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
	}
    }
    #cat("\n")#- not in S
    invisible(x)
}

###Do the F tests for a cablibrated GWR model
F1234.test<-function(F.test.parameters=list())
{
  F1.test<-F2.test<-F3.test<-F4.test<-NULL
	v1 <-F.test.parameters$tr.S
  v2 <- F.test.parameters$tr.StS
  dp.n <- F.test.parameters$dp.n
  #effective d.f. is n - 2*v1 + v2
	edf <- dp.n - 2*v1 + v2
  #####################################F1,F2,F3 test according to Leung et al. 2000
  RSSg <- F.test.parameters$RSS.gw
  RSSo <- F.test.parameters$RSS.lm
  DFo <- F.test.parameters$DF.lm
  ######distance matrix calculaion bug by Jasdeep
  dMat <- F.test.parameters$dMat
  if(dim(dMat)[1]==1)
  {
     DM.given <- F
     dp.locat <- F.test.parameters$dp.locat
     p <- F.test.parameters$p
     theta <- F.test.parameters$theta
     longlat <- F.test.parameters$longlat
  }
  else
    DM.given <- T
  delta1 <- dp.n-2*v1+v2
  q.diag <- F.test.parameters$q.diag
  sigma2.delta1 <- RSSg/delta1
  sigma2 <- RSSg/dp.n
  odelta2 <- sum(q.diag^2)
  var.n<- F.test.parameters$var.n
  x <- F.test.parameters$x
  betas <- F.test.parameters$betas
  bw <- F.test.parameters$bw
  adaptive <- F.test.parameters$adaptive
  kernel <- F.test.parameters$kernel
 	delta2 <- 0
  L.delta1 <- sum(q.diag)
  L.delta2 <- delta2
  #####F1 test
  F1<-(RSSg/L.delta1)/(RSSo/DFo)
  F1.DF<-c(L.delta1^2/L.delta2, DFo)
  F1.pv <- pf(F1, F1.DF[1], F1.DF[2], lower.tail=TRUE)
  F1.test<-matrix (nrow = 1, ncol = 4)
  F1.test[1,1]<-F1
  F1.test[1,c(2,3)]<-F1.DF
  F1.test[1,4]<-F1.pv
  colnames(F1.test) <- c(
    "F1 statistic",
    "Numerator DF",
    "Denominator DF",
    "Pr(>)"
  )
  rownames(F1.test)<-c(" ")
  #####F2 Test
  F2<-((RSSo-RSSg)/(DFo-L.delta1))/(RSSo/DFo)
  F2.DF<-c((DFo-L.delta1)^2/(DFo-2*L.delta1+L.delta2),DFo)
  F2.pv <- pf(F2, F2.DF[1], F2.DF[2], lower.tail=FALSE)
  F2.test<-matrix (nrow = 1, ncol = 4)
  F2.test[1,1]<-F2
  F2.test[1,c(2,3)]<-F2.DF
  F2.test[1,4]<-F2.pv
  colnames(F2.test) <- c(
    "F2 statistic",
    "Numerator DF",
    "Denominator DF",
    "Pr(>)"
  )
  rownames(F2.test)<-NULL
  ####F3 Test
  ek <- diag (var.n)
  iden <- diag (dp.n)
  J <- matrix (1, nrow = dp.n, ncol = dp.n)
  ###Vk2: Vk squre
  Vk2 <- numeric (var.n)
  for (i in 1:var.n)
  {
	   Vk2[i] <- (1/dp.n)*(t(betas[,i])%*%(iden-(1/dp.n)*J)%*%betas[,i])
  }
  gamma1 <- numeric (var.n)
  gamma2 <- numeric (var.n)
  numdf <- numeric (var.n)
  F3 <- numeric (var.n)
  F3.pv<- numeric (var.n)
  F3.DF<-matrix(numeric (var.n*2),ncol=2)
  for (i in 1:var.n)
  {
	   B <- matrix (nrow = dp.n, ncol = dp.n)
	   for (j in 1:dp.n)
      {
	     if(DM.given)
	     {
	       dist.vj<-dMat[,j]
	     }
	     else
	     {
	       dist.vj<- gw.dist(dp.locat,dp.locat, focus=j, p, theta, longlat)
	     }
        
		    wj <- gw.weight(dist.vj,bw,kernel,adaptive)
		    B[j,] <- ek[i,] %*% solve(t(x)%*%diag(wj)%*%x) %*%t(x) %*% diag(wj)
	   }
	   BJ<- (1/dp.n)*(t(B)%*%(iden-(1/dp.n)*J)%*%B)
	   gamma1[i] <- sum(diag(BJ))
   	 gamma2[i] <- sum(diag(BJ)^2)
	   numdf[i] <- gamma1[i]^2/gamma2[i]
	   #F3[i] <- (Vk2[i]/gamma1[i])/sigma2
     F3[i] <- (Vk2[i]/gamma1[i])/sigma2.delta1
	   F3.pv[i] <- pf(F3[i], numdf[i], F1.DF[1], lower.tail=FALSE)
   F3.DF[i,]<-c(numdf[i], F1.DF[1])
  }
  F3.test <- matrix (nrow = var.n, ncol = 4)
  F3.test[,1] <- F3
  F3.test[,c(2,3)] <- F3.DF
	F3.test[,4] <- F3.pv
  colnames(F3.test) <- c(
    "F3 statistic",
    "Numerator DF",
    "Denominator DF",
    "Pr(>)"
  )
  rownames(F3.test) <- colnames(x)
  ######################F4, the test statistic adopted f rom the GWR book p92
  F4<-RSSg/RSSo
  F4.DF<-c(delta1, DFo)
  F4.pv<-pf(F4, F4.DF[1], F4.DF[2], lower.tail=TRUE)
  F4.test<-matrix (nrow = 1, ncol = 4)
  F4.test[1,1]<-F4
  F4.test[1,c(2,3)]<-F4.DF
  F4.test[1,4]<-F4.pv
  colnames(F4.test) <- c(
    "F4 statistic",
    "Numerator DF",
    "Denominator DF",
    "Pr(>)"
  )
  rownames(F4.test)<-NULL
  res<-list(F1.test=F1.test,F2.test=F2.test,F3.test=F3.test,F4.test=F4.test)
  res
}

# Randomisation tests for GWR parameters

test.gwr.par<-function(formula, data, regression.points, bw, kernel = "bisquare",
    adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat,W.vect=NULL,nperm=999)
{
	timings <- list()
    timings[["start"]] <- Sys.time()
    this.call <- match.call()
    p4s <- as.character(NA)
    if (missing(regression.points)) {
        rp.given <- FALSE
        regression.points <- data
        hatmatrix <- T
    }
    else {
        rp.given <- TRUE
        hatmatrix <- F
    }
    if (is(data, "Spatial")) {
        p4s <- proj4string(data)
        dp.locat <- coordinates(data)
        data <- as(data, "data.frame")
    }
    else {
        stop("Given regression data must be Spatial*DataFrame")
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    var.n <- ncol(x)
    rp.locat <- coordinates(regression.points)
    rp.n <- nrow(rp.locat)
    dp.n <- nrow(data)
    betas <- matrix(nrow = rp.n, ncol = var.n)
    betas.SE <- matrix(nrow = rp.n, ncol = var.n)
    betas.TV <- matrix(nrow = rp.n, ncol = var.n)
    S <- matrix(nrow = dp.n, ncol = dp.n)
    colnames(betas) <- colnames(x)
    colnames(betas)[1] <- "Intercept"
    lms <- lm(formula, data = data)
    lms$x <- x
    lms$y <- y
    gTSS <- c(cov.wt(matrix(y, ncol = 1), wt = rep(as.numeric(1),
        dp.n), method = "ML")$cov * dp.n)
    if (missing(dMat))
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != rp.n)
            stop("Dimensions of dMat are not correct")
    }

    beta.i<-array(dim=c(dp.n,2,nperm+1))

    for (i in 1:rp.n) {
        if (DM.given)
            dist.vi <- dMat[, i]
        else {
            if (rp.given) dist.vi <- gw.dist(dp.locat, rp.locat, focus = i, p, theta, longlat)
            else dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, p = p, theta = theta, longlat = longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)

 if(!is.null(W.vect)) W.i<-W.i*W.vect

## Original Data

        gwsi <- gw_reg(x, y, W.i, hatmatrix, i)
        bbeta<-gwsi[[1]]
        beta.i[i, ,1]<-bbeta
###
	for(perm in 1:nperm)
	{
		a<-sample(seq(1,dp.n))
		xx<-x[a,]
		yy<-y[a]

        gwsi <- gw_reg(xx, yy, W.i, hatmatrix, i)
        beta.i[i, ,perm+1] <- gwsi[[1]]
    }
  }

#############

	dimnames(beta.i)[[1]]<-seq(0,dp.n-1)
	dimnames(beta.i)[[2]]<-rownames(bbeta)
	dimnames(beta.i)[[3]]<-c('Original_Data',paste('Sample',seq(1:nperm),sep='_'))

	test<-apply(beta.i,c(1,2),function(x,n) rank(x,ties.method='first')[1]/n,n=nperm+1)

	test
}
