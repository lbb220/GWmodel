####Scalable GWR
gwr.scalable <- function(formula, data, bw.adapt=100, kernel = "gaussian", polynomial = 4, p = 2, theta = 0, longlat = F, dMat) 
{
  timings <- list()
  timings[["start"]] <- Sys.time()
  # # Match Variables
  this.call <- match.call()
  p4s <- as.character(NA)
  # # Check the given data frame and regression points
  hatmatrix<-T
  # ## Data points
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  } else {
    stop("Given regression data must be Spatial*DataFrame")
  }
  
  # # Extract the data frame
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
  #dp.n <- nrow(data)
  dp.n <- nrow(x)
  rp.n <- dp.n
  dmat.p <- p
  betas <- matrix(0, rp.n, var.n)
  if (hatmatrix) {
    betas.SE <- matrix(0, rp.n, var.n)
    betas.TV <- matrix(0, rp.n, var.n)
  }
  idx1 <- match("(Intercept)", colnames(x))
  if (!is.na(idx1)) colnames(x)[idx1] <- "Intercept"
  colnames(betas) <- colnames(x)
  xnames <- colnames(x)
  
  # # Linear regression
  lms <- lm(formula = formula, data = data)
  lms$x <- x
  lms$y <- y
  gTss <- c(cov.wt(matrix(y, ncol = 1), wt = rep(as.numeric(1), dp.n), method = "ML")$cov * dp.n)
  
  # # Distance Mat is given or not
  neighbours <- numeric(bw.adapt^2)
  neighbours.dist <- NULL
  DM.given <- F
  DM1.given <- F
  if (missing(dMat)) {
    DM.given <- F
    DM1.given <- F
    if (dmat.p == 2 & theta == 0 & longlat == F) {
      data.knn <- get.knn(dp.locat, k = bw.adapt)
      neighbours <- data.knn$nn.index
      neighbours.dist <- data.knn$nn.dist
      rm(data.knn)
      gc()
    } else {
      if (dp.n <= 10000 & dmat.p == 2 & theta == 0) {
        neighbours <- knearneigh(dp.locat, k = bw.adapt, longlat = longlat)$nn
        dMat <- gw.dist(dp.locat = dp.locat, p = 2, theta = theta, longlat = longlat)
        DM.given <- T
      } else {
        neighbours <- matrix(0, dp.n, bw.adapt)
        neighbours.dist <- matrix(0, dp.n, bw.adapt)
        time0 <- Sys.time()
        for (i in 1:dp.n) {
          dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, p = dmat.p, theta = theta, longlat = longlat)
          dist.vi.sort <- sort(dist.vi, index.return = TRUE)
          neighbours[i,] <- dist.vi.sort$ix[2:(bw.adapt + 1)]
          neighbours.dist[i,] <- dist.vi.sort$x[2:(bw.adapt + 1)]
        }
      }
    }
  } else {
    DM.given <- T
    DM1.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != rp.n) {
      stop("Dimensions of dMat are not correct!")
    }
  }
  if (DM.given) {
    dMat.sort <- apply(dMat, 2, sort, index.return = TRUE)
    neighbours <- t(sapply(dMat.sort, function (i) {i$ix}, simplify = "array"))[,2:(bw.adapt + 1)]
    neighbours.dist <- t(sapply(dMat.sort, function (i) {i$x}, simplify = "array"))[,2:(bw.adapt + 1)]
  }
  timings[["precalibration"]] <- Sys.time()
  # # Calibration the model
  P <- polynomial
  band0 <- 0
  G0 <- 1
  # ## Base bandwidth and base kernel
  if (tolower(kernel) == "gaussian") {
    band0 <- median(neighbours.dist[, min(50, bw.adapt)]) / sqrt(3)
    G0 <- exp(-(neighbours.dist / band0)^2)
  } else if (tolower(kernel) == "exponential") {
    band0 <- median(neighbours.dist[, min(50, bw.adapt)]) / 3
    G0 <- exp(-(neighbours.dist / band0)^2)
  }
  XtX <- crossprod(x)
  XtY <- crossprod(x, y)
  
  M0 <- scgwr_pre(x, y, bw.adapt, P, band0, G0, neighbours)
  Mx0 <- M0$Mx0
  My0 <- M0$My0
  timings[["precalculate"]] <- Sys.time()
  
  b.tilde <- 1
  alpha <- 0.01
  optim.result <- optim(par = c(b.tilde, alpha), fn = scgwr_loocv, gr = NULL, x, y, bw.adapt, P, Mx0, My0, XtX, XtY)
  b.tilde <- optim.result$par[1]^2
  alpha <- optim.result$par[2]^2
  prameters <- c(b.tilde, alpha)
  timings[["optimize"]] <- Sys.time()
  
  cv.score <- optim.result$value
  result.reg <- scgwr_reg(x, y, bw.adapt, P, G0, Mx0, My0, XtX, XtY, neighbours, prameters)
  betas <- result.reg$betas
  betas.SE <- result.reg$betas.SE
  tr.S <- result.reg$tr.S
  tr.StS <- result.reg$tr.StS
  timings[["calibration"]] <- Sys.time()
  # # Diagnostic statistics
  yhat <- rowSums(x * betas)
  residual <- y - yhat
  RSS <- c(sum(residual^2))
  
  enp <- 2 * tr.S - tr.StS
  edf <- dp.n - 2 * tr.S + tr.StS
  sig <- c(sqrt(RSS / (dp.n - enp)))
  betas.SE <- sig * betas.SE
  bt <- betas / betas.SE
  bp0 <- data.frame(2 - 2 * pt(abs(bt), df = c(dp.n) - enp))
  bp <- data.frame(bp0 * c(abs(enp / var.n)))
  bp[bp > 1] <- 1
  betas <- data.frame(betas)
  betas.SE <- data.frame(betas.SE)
  betas.TV <- data.frame(betas / betas.SE)
  names(betas) <- names(betas.SE) <- names(betas.TV) <- names(bp0) <- names(bp) <- xnames
  spar <- data.frame(c(b.tilde, alpha))
  names(spar) <- "Estimate"
  rownames(spar) <- c("scale", "penalty")
  loglik <- -dp.n * log(sig) - dp.n / 2 * log(2 * pi)
  AIC <- -2 * loglik + dp.n + tr.S
  AICc <- -2 * loglik + dp.n * (dp.n + tr.S) / (dp.n - 2 - tr.S)
  R2 <- cor(yhat, y)^2
  adjR2 <- 1 - (1 - R2) * (dp.n - 1)/(dp.n - enp - 1)
  CV <- sqrt(cv.score / dp.n)
  
  # # Generate return result
  # ## GWR arguments
  GW.arguments <- list(
    formula = formula,
    data = data, 
    bw = bw.adapt, 
    kernel = kernel, 
    polynomial = polynomial, 
    p = dmat.p, 
    theta = theta, 
    longlat = longlat,
    DM.given = DM.given,
    hatmatrix = TRUE,
    F123.test = FALSE
  )
  # ## GWR diagnostic
  GW.diagnostic <- list(
    AIC = AIC,
    AICc = AICc,
    enp = enp,
    edf = edf,
    gw.R2 = R2,
    gwR2.adj = adjR2,
    RSS.gw = RSS
  )
  # ## GWR result Data Frame
  gwres.df <- data.frame(betas, y, yhat, residual, CV, betas.SE, betas.TV)
  colnames(gwres.df) <- c(c(c(colnames(betas), c("y", "yhat", "residual", "CV_Score")),
                            paste(colnames(betas), "SE", sep = "_")),
                          paste(colnames(betas), "TV", sep = "_"))
  # ## SDF
  if (is(data, "SpatialPolygonsDataFrame")) {
    polygons <- polygons(data)
    rownames(gwres.df) <- sapply(slot(polygons, "polygons"), function (i) slot(i, "ID"))
    SDF <- SpatialPolygonsDataFrame(Sr = polygons, data = gwres.df, match.ID = F)
  } else {
    SDF <- SpatialPointsDataFrame(coords=dp.locat, data=gwres.df, proj4string = CRS(p4s), match.ID = F)
  }
  timings[["stop"]] <- Sys.time()
  res <- list(
    GW.arguments = GW.arguments,
    GW.diagnostic = GW.diagnostic,
    lm = lms,
    SDF = SDF,
    this.call = this.call,
    timings = timings
  )
  class(res) <- "scgwrm"
  invisible(res)
}
gwr.scalable.loocv <- function(target, params) {
  # # Get parameters
  pos.px <- params$pos.px
  pos.var1 <- params$pos.var1
  pos.var2 <- params$pos.var2
  pos.py <- params$pos.py
  pos.vary <- params$pos.vary
  Mx0 <- params$Mx0
  My0 <- params$My0
  P <- params$P
  var.n <- params$var.n
  dp.n <- params$dp.n
  x <- params$x
  y <- params$y
  XtX <- params$XtX
  XtY <- params$XtY
  # # optimization
  b <- target[1]^2
  alpha <- target[2]^2
  R <- rep(b, P + 1)
  for (p in 2:(P + 1)) {
    R[p] <- b^p
  }
  Mx <- R[pos.px] * Mx0
  My <- R[pos.py] * My0
  Mx.sum <- matrix(0, var.n * var.n, dp.n)
  My.sum <- matrix(0, var.n, dp.n)
  for (k1 in 1:var.n) {
    for (k2 in 1:var.n) {
      Mx.sum[(k1 - 1)*var.n + k2, ] <- colSums(Mx[(pos.var1 == k1) & (pos.var2 == k2), ])
    }
    My.sum[k1, ] <- colSums(My[(pos.vary == k1), ])
  }
  Mx.list <- list(NULL)
  My.list <- list(NULL)
  for (i in 1:dp.n) {
    Mx.list[[i]] <- matrix(Mx.sum[, i], var.n, var.n) + alpha * XtX
    My.list[[i]] <- My.sum[, i] + alpha * XtY
  }
  tryres <- try(Mx.inv <- lapply(Mx.list, solve), silent = TRUE)
  if (class(tryres) == "try-error") {
    error <- 10^10
  } else {
    beta <- t(matrix(unlist(mapply("%*%", Mx.inv, My.list, SIMPLIFY = FALSE)), nrow = var.n, ncol = dp.n))
    yhat <- rowSums(x * beta)
    error <- sum((y - yhat)^2)
  }
  return(error)
}

print.scgwrm <- function(x, ...)
{
  if(class(x) != "scgwrm") stop("It's not a scalable gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars <- all.vars(x$GW.arguments$formula)
  var.n <- length(x$lm$coefficients)
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
  cat("\n   ***********************************************************************\n")
  cat("   Program stops at:", as.character(x$timings$stop), "\n")
  invisible(x)
}