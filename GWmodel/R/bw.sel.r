###Basic function for bandwidth selction
##Author: Binbin Lu
bw.gwr<-function(formula, data, approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, 
                longlat=F,dMat,parallel.method=F,parallel.arg=NULL)
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
  ####################################Joke
  if(dp.n>1500)
  {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("          -----A kind suggestion from GWmodel development group\n")
   }  
  #################### Recommond to specify a distance matrix
  if (missing(dMat))
  {
      dMat <- NULL
      DM.given<-F
      if(dp.n + dp.n <= 10000)
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
    stop ("Dimensions of dMat are not correct")
  }
  #########Find the range of the fixed bandwidth
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
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
  bw<-NA
  # ## make cluseter
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) {
      cl.n <- max(detectCores() - 4, 2)
      parallel.arg <- makeCluster(cl.n)
    } else cl.n <- length(parallel.arg)
    clusterCall(parallel.arg, function() { library(GWmodel) })
  }
  # ## call for functions
  if(approach == "bic" || approach == "BIC")
      bw <- gold(gwr.bic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p, theta, longlat, dMat, T, parallel.method, parallel.arg)
  else if(approach == "aic" || approach == "AIC" || approach == "AICc")
      bw <- gold(gwr.aic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p, theta, longlat, dMat, T, parallel.method, parallel.arg)    
  else 
      bw <- gold(gwr.cv, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p, theta, longlat, dMat, T, parallel.method, parallel.arg)
  # ## stop cluster
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) stopCluster(parallel.arg)
  }
  bw

}





####Calculate the CV score with a given bandwidth
##Author: Binbin Lu
gwr.cv<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T,parallel.method=F,parallel.arg=NULL)
{
  dp.n<-length(dp.locat[,1])
  #########Distance matrix is given or not

  if (missing(dMat)) dMat <- matrix(0, 0, 0)
  else if (is.null(dMat) || !is.matrix(dMat)) {
    DM.given<-F
    dMat <- matrix(0, 0, 0)
  }
  else {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################CV
  gw.resi <- NULL
  if (parallel.method == FALSE) {
    gw.resi <- try(gw_cv_all(X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive))
    if(!inherits(gw.resi, "try-error")) CV.score <- gw.resi 
    else CV.score <- Inf
  } else if (parallel.method == "omp") {
    if (missing(parallel.arg)) { threads <- 0 } else {
      threads <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    gw.resi <- try(gw_cv_all_omp(X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive, threads))
    if(!inherits(gw.resi, "try-error")) CV.score <- gw.resi 
    else CV.score <- Inf
  } else if (parallel.method == "cuda") {
    if (missing(parallel.arg)) { groupl <- 0 } else {
      groupl <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    gw.resi <- try(gw_cv_all_cuda(X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive, groupl))
    if(!inherits(gw.resi, "try-error")) CV.score <- gw.resi 
    else CV.score <- Inf
  } else if (parallel.method == "cluster") {
    print("Parallel using cluster.")
    cl.n <- length(parallel.arg)
    cl.results <- clusterApplyLB(parallel.arg, 1:cl.n, function(group.i, cl.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive) {
      cv.result <- try(gw_cv_all(x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive, cl.n, group.i))
      if(!inherits(gw.resi, "try-error")) return(cv.result) 
      else return(Inf)
    }, cl.n, X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive)
    gw.resi <- unlist(cl.results)
    if (!all(is.infinite(gw.resi))) CV.score <- sum(gw.resi)
    else CV.score <- Inf
  } else {
    CV<-numeric(dp.n)
    for (i in 1:dp.n) {
      if (DM.given) dist.vi<-dMat[,i]
      else dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      W.i[i]<-0
      gw.resi<- try(gw_reg(X, Y, W.i, FALSE, i))
      if(!inherits(gw.resi, "try-error")) {
        yhat.noi<-X[i,]%*%gw.resi[[1]]
        CV[i]<-Y[i]-yhat.noi
      } else {
        CV[i]<-Inf
        break
      }
    }
    if (!any(is.infinite(CV))) CV.score<-t(CV) %*% CV
    else CV.score<-Inf
  }
  if(verbose) {
    if(adaptive) cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
    else cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
  }
  ifelse(is.nan(CV.score), Inf, CV.score)
}

gwr.cv.contrib<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat,parallel.method=F,parallel.arg=NULL)
{
   dp.n<-length(dp.locat[,1])
   #########Distance matrix is given or not

  if (is.null(dMat)) DM.given<-F
  else {
    if(dim(dMat)[1]==1) {
      DM.given <- F
    } else {
      DM.given<-T
      dim.dMat<-dim(dMat)
      if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
        stop("Dimensions of dMat are not correct")
    }
  }
  ############################################CV
  CV<-numeric(dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    #W.i<-gwr.Gauss(dist.vi^2, bw)
    #print(W.i)
    W.i[i]<-0
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    #fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    #gw.resi<-try(fun1(X,Y,W.i))
    gw.resi<- try(gw_reg(X, Y, W.i, FALSE, i))
  
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))

    if(!inherits(gw.resi, "try-error"))
    {
      yhat.noi <- X[i,] %*% gw.resi[[1]]
      CV[i] <- Y[i] - yhat.noi
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  CV
}
####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
gwr.aic<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T,parallel.method=F,parallel.arg=NULL)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  #########Distance matrix is given or not

  if (missing(dMat)) {
    DM.given <- F
    dMat <- matrix(0, 0, 0)
  }
  else if (is.null(dMat) || !is.matrix(dMat)) {
    DM.given<-F
    dMat <- matrix(0, 0, 0)
  }
  else {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  AICc.value <- Inf
  if (parallel.method == FALSE) {
    res <- try(gw_reg_all(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      AICc.value <- AICc1(Y, X, betas, s_hat)
    }
    else AICc.value <- Inf
  } else if (parallel.method == "omp") {
    if (missing(parallel.arg)) { threads <- 0 } else {
      threads <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    res <- try(gw_reg_all_omp(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, threads))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      AICc.value <- AICc1(Y, X, betas, s_hat)
    }
    else AICc.value <- Inf
  } else if (parallel.method == "cuda") {
    if (missing(parallel.arg)) { groupl <- 0 } else {
      groupl <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    res <- try(gw_reg_all_cuda(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, groupl))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      AICc.value <- AICc1(Y, X, betas, s_hat)
    }
    else AICc.value <- Inf
  } else if (parallel.method == "cluster") {
    cl.n <- length(parallel.arg)
    cl.results <- clusterApplyLB(parallel.arg, 1:cl.n, function(group.i, cl.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive) {
      res <- try(gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, cl.n, group.i))
      if(!inherits(res, "try-error")) return(res) 
      else return(NULL)
    }, cl.n, X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive)
    if (!any(is.null(cl.results))) {
      betas <- matrix(0, nrow = dp.n, ncol=var.n)
      s_hat <- numeric(2)
      for (i in 1:cl.n) {
        res <- cl.results[[i]]
        betas = betas + res$betas
        s_hat = s_hat + res$s_hat
      }
      AICc.value <- AICc1(Y, X, betas, s_hat)
    } else AICc.value <- Inf
  } else {
    s_hat <- numeric(2)
    betas <- matrix(nrow = dp.n, ncol = var.n)
    for (i in 1:dp.n) {
      if (DM.given) dist.vi <- dMat[,i]
      else dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
      W.i <- gw.weight(dist.vi,bw,kernel,adaptive)
      res<- try(gw_reg(X,Y,W.i,TRUE,i))
      if(!inherits(res, "try-error")) {
        si <- res[[2]]
        s_hat[1] = s_hat[1] + si[i]
        s_hat[2] = s_hat[2] + sum(tcrossprod(si))
        betas[i,] <- res[[1]]
      } else {
        s_hat[1] <- Inf
        s_hat[2] <- Inf
        break
      }  
    }
    if (!any(is.infinite(s_hat))) {
      AICc.value <- AICc1(Y, X, betas, s_hat)
    } else AICc.value<-Inf
  }
  if(is.nan(AICc.value)) AICc.value <- Inf
  if(verbose) {     
    if(adaptive) cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc.value, "\n")
    else cat("Fixed bandwidth:", bw, "AICc value:", AICc.value, "\n")
  }
  AICc.value
}

####Calculate the BIC with a given bandwidth
##Author: Binbin Lu
gwr.bic<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T,parallel.method=F,parallel.arg=NULL)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  #########Distance matrix is given or not

  if (missing(dMat)) {
    DM.given <- F
    dMat <- matrix(0, 0, 0)
  }
  else if (is.null(dMat) || !is.matrix(dMat)) {
    DM.given<-F
    dMat <- matrix(0, 0, 0)
  }
  else {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  BIC.value <- Inf
  if (parallel.method == FALSE) {
    res <- try(gw_reg_all(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      BIC.value <- BIC(Y, X, betas, s_hat)
    }
    else BIC.value <- Inf
  } else if (parallel.method == "omp") {
    if (missing(parallel.arg)) { threads <- 0 } else {
      threads <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    res <- try(gw_reg_all_omp(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, threads))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      BIC.value <- BIC(Y, X, betas, s_hat)
    }
    else BIC.value <- Inf
  } else if (parallel.method == "cuda") {
    if (missing(parallel.arg)) { groupl <- 0 } else {
      groupl <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
    }
    res <- try(gw_reg_all_cuda(X, Y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, groupl))
    if(!inherits(res, "try-error")) {
      betas <- res$betas
      s_hat <- res$s_hat
      BIC.value <- BIC(Y, X, betas, s_hat)
    }
    else BIC.value <- Inf
  } else if (parallel.method == "cluster") {
    cl.n <- length(parallel.arg)
    cl.results <- clusterApplyLB(parallel.arg, 1:cl.n, function(group.i, cl.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive) {
      res <- try(gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, cl.n, group.i))
      if(!inherits(res, "try-error")) return(res) 
      else return(NULL)
    }, cl.n, X, Y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive)
    if (!any(is.null(cl.results))) {
      betas <- matrix(0, nrow = dp.n, ncol=var.n)
      s_hat <- numeric(2)
      for (i in 1:cl.n) {
        res <- cl.results[[i]]
        betas = betas + res$betas
        s_hat = s_hat + res$s_hat
      }
      BIC.value <- BIC(Y, X, betas, s_hat)
    } else BIC.value <- Inf
  } else {
    s_hat <- numeric(2)
    betas <- matrix(nrow = dp.n, ncol = var.n)
    for (i in 1:dp.n) {
      if (DM.given) dist.vi <- dMat[,i]
      else dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
      W.i <- gw.weight(dist.vi,bw,kernel,adaptive)
      res<- try(gw_reg(X,Y,W.i,TRUE,i))
      if(!inherits(res, "try-error")) {
        si <- res[[2]]
        s_hat[1] = s_hat[1] + si[i]
        s_hat[2] = s_hat[2] + sum(tcrossprod(si))
        betas[i,] <- res[[1]]
      } else {
        s_hat[1] <- Inf
        s_hat[2] <- Inf
        break
      }  
    }
    if (!any(is.infinite(s_hat))) {
      BIC.value <- BIC(Y, X, betas, s_hat)
    } else BIC.value <-Inf
  }
  if(is.nan(BIC.value)) BIC.value <- Inf
  if(verbose) {     
    if(adaptive) cat("Adaptive bandwidth (number of nearest neighbours):", bw, "BIC value:", BIC.value, "\n")
    else cat("Fixed bandwidth:", bw, "BIC value:", BIC.value, "\n")
  }
  BIC.value
}
#######################################################################
#################### Golden Section search ############################
#######################################################################
# Golden selection function 
###Author: created by MC, edited by BL
# from www.me.gatech.edu/~aferri/me2016/Golden_Section.ppt
gold<-function(fun,xL,xU,adapt.bw=F,...){
   eps=1e-4 # working value - can be changed
   R <- (sqrt(5)-1)/2  # R <- 0.61803398....
   iter <- 1
   d <- R*(xU-xL)
   if (adapt.bw)
   {
     x1 <- floor(xL+d)
     x2 <- round(xU-d)
   }
   else
   {
     x1 <- xL+d
     x2 <- xU-d
   }  
   f1 <- eval(fun(x1,...))
   f2 <- eval(fun(x2,...))
   d1<-f2-f1
# Establish initial value of xopt:
   if (f1 < f2)
      xopt  <-  x1
   else xopt  <-  x2
# start main loop
   ea <- 100
   while ((abs(d) > eps) && (abs(d1) > eps)) {
      d <- R*d
      if   (f1 < f2) {
         xL  <-  x2
         x2  <-  x1
         if (adapt.bw)         
           x1 <- round(xL+d)
         else
           x1 <- xL+d
         #x1  <-  xL + d
         f2  <-  f1
         f1 <- eval(fun(x1,...))
         }
      else {
         xU  <-  x1
         x1  <-  x2
         if (adapt.bw)         
           x2 <- floor(xU - d)
         else
           x2  <-  xU - d 
         f1  <-  f2
         f2 <- eval(fun(x2,...))
         }
   iter  <-  iter + 1
# Establish value of xopt after iteration:
   if    (f1 < f2)
     xopt  <-  x1
   else xopt  <-  x2
   d1<-f2-f1
   ##print(paste(iter,f1,f2,xopt))
   }
    xopt
}
