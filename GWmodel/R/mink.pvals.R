#################################################################
#Select values of p automatically
gwr.mink.pval <- function(formula, data, criterion="AIC", bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", left.interval=0.25,
                       right.interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given regression data must be Spatial*DataFrame")
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  var.n <- ncol(x)
  longlat <- F
  bw.given <- T
  if(missing(bw))
     bw.given <- F
  #############
  p <- 2
  p.vals <- c(p)
  cretion.vals <- c()
  cat("  ==========  Try MinkoWski distance with p: ", p)
  cat("  ==========\n")
  dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta0, longlat=longlat)
  if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
  if(criterion=="AIC" || criterion=="AICc")
  {
    res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
    cretion.eu <- res1[[1]]
    cretion.vals <- c(cretion.eu)
    if(adaptive)
      cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", cretion.eu, "\n")
    else
      cat("    Fixed bandwidth: ",bw, ", AICc value: ", cretion.eu, "\n")
  }
  else
  {
    res1 <- gwr.cv1(bw, x, y, kernel,adaptive, dp.locat, dMat)
    cretion.eu <- res1[[1]]
    cretion.vals <- c(cretion.eu)
    if(adaptive)
      cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", CV score: ", cretion.eu, "\n")
    else
      cat("    Fixed bandwidth: ",bw, ", CV score: ", cretion.eu, "\n")
  }
  #####################Larger than 2 by the right.interval
  #Calculate the cretion value when p=inf
  cat("  ==========  Try Minkovski distance with p: ", Inf)
  cat("  ==========\n")
  dMat <- gw.dist(dp.locat=dp.locat, p=Inf, theta=theta0, longlat=longlat)
  if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
  if(criterion=="AIC" || criterion=="AICc")
  {
    res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
    cretion.inf <- res1[[1]]
    if(adaptive)
      cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", cretion.inf, "\n")
    else
      cat("    Fixed bandwidth: ",bw, ", AICc value: ", cretion.inf, "\n")
  }
  else
  {
    res1 <- gwr.cv1(bw, x, y, kernel,adaptive, dp.locat, dMat)
    cretion.inf <- res1[[1]]
    if(adaptive)
      cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", CV score: ", cretion.inf, "\n")
    else
      cat("    Fixed bandwidth: ",bw, ", CV score: ", cretion.inf, "\n")
  } 
  #######Start the search from 2
  cretion.val <- cretion.eu
  p.right <- 2 + right.interval 
  while(abs(cretion.val-cretion.inf)>drop.tol)
  {
    p.vals <- c(p.vals, p.right)
    cat("  ==========  Try MinkoWski distance with p: ", p.right)
    cat("  ==========\n")
    dMat <- gw.dist(dp.locat=dp.locat, p=p.right, theta=theta0, longlat=longlat)
    print(dMat[1,2])
    print(p.right)
    if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
    if(criterion=="AIC" || criterion=="AICc")
    {
      res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
      cretion.val <- res1[[1]]
      cretion.vals <- c(cretion.vals, cretion.val)
      if(adaptive)
        cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", cretion.val, "\n")
      else
        cat("    Fixed bandwidth: ",bw, ", AICc value: ", cretion.val, "\n")
    }
    else
    {
      res1 <- gwr.cv1(bw, x, y, kernel,adaptive, dp.locat, dMat)
      cretion.val <- res1[[1]]
      cretion.vals <- c(cretion.vals, cretion.val)
      if(adaptive)
        cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", CV score: ", cretion.val, "\n")
      else
        cat("    Fixed bandwidth: ",bw, ", CV score: ", cretion.val, "\n")
    }
    p.right <- p.right + right.interval   
  }
  #Less than 2 by the left.interval
  p.left <- 2 - left.interval
  while(p.left>0)
  {
    p.vals <- c(p.left, p.vals)
    cat("  ==========  Try Minkovski distance with p: ", p.left)
    cat("  ==========\n")
    dMat <- gw.dist(dp.locat=dp.locat, p=p.left, theta=theta0, longlat=longlat)
    if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
    if(criterion=="AIC" || criterion=="AICc")
    {
      res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
      cretion.val <- res1[[1]]
      cretion.vals <- c(cretion.val, cretion.vals)
      if(adaptive)
        cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", cretion.val, "\n")
      else
        cat("    Fixed bandwidth: ",bw, ", AICc value: ", cretion.val, "\n")
    }
    else
    {
      res1 <- gwr.cv1(bw, x, y, kernel,adaptive, dp.locat, dMat)
      cretion.val <- res1[[1]]
      cretion.vals <- c(cretion.val, cretion.vals)
      if(adaptive)
        cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", CV score: ", cretion.val, "\n")
      else
        cat("    Fixed bandwidth: ",bw, ", CV score: ", cretion.val, "\n")
    }
    p.left <- p.left - left.interval   
  }
  p.vals <- c(p.vals, Inf)
  cretion.vals <- c(cretion.vals, cretion.inf)
  n <- length(p.vals)
  p.dropped <- c(FALSE)
   cretion.val0 <- cretion.vals[1]
   for(i in 2:n)
   {
     cretion.val1 <- cretion.vals[i]
     if(p.vals[i]==2)
     {
       p.dropped <- c(p.dropped ,FALSE)
       cretion.val0 <- cretion.val1
       next
     }
     if(abs(cretion.val1-cretion.val0)>drop.tol)
     {
        p.dropped <- c(p.dropped ,FALSE)
        cretion.val0 <- cretion.val1
     }
     else
        p.dropped <- c(p.dropped ,TRUE)
   }
   #vis.df <- data.frame(as.factor(p.vals), cretion.vals, p.dropped)
   #rownames(vis.df) <- as.character(p.vals)
   #colnames(vis.df) <- c("pvals", "AICc", "Dropped")
   #dotplot(p.vals~AICc, vis.df)
  res <- list(p.vals=p.vals, cretion.vals=cretion.vals, p.dropped=p.dropped)
  class(res) <- "pvlas"
  res
}
plot.pvlas <- function(x, ...)
{
   p.vals <- x[[1]]
   cretion.vals <- x[[2]]
   p.dropped <- x[[3]] 
   x.offset <- 0.5
   y.offset <- 1
   n <- length(p.vals)
   idx <- 1:n
   plot(idx, cretion.vals, typ="l", lty=2, xaxt='n', xlab="value of p", ylab="AICc value")
   points(idx[!p.dropped], cretion.vals[!p.dropped], pch=16,col="blue")
   points(idx[p.dropped], cretion.vals[p.dropped], pch=4,col="red")
   text(idx[!p.dropped]+x.offset,cretion.vals[!p.dropped]+y.offset, p.vals[!p.dropped], font=2,cex=0.75)
   #text(idx[p.dropped]+0.5,cretion.vals[p.dropped]+1, p.vals[p.dropped], font=2,col="grey",cex=0.75)
   legend("topright",legend=c("Values of p to be used", "Values of p to be dropped"),col=c("blue", "red"), pch=c(16,4),)
    
}


gwr.mink.pval.forward <- function(formula, data, bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", p.max=Inf,p.min=2,
                       interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given regression data must be Spatial*DataFrame")
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  var.n <- ncol(x)
  longlat <- F
  bw.given <- T
  if(missing(bw))
     bw.given <- F
   
  #############
  # the AICc value calculated from the maximum value of p
  AIC.pmax <- 0
  p.vals <- c()
  AICc.vals <- c()
  p.dropped <- c()
  bws <- c()
  
  cat("  ==========  Try MinkoWski distance with the maximum p: ", p.max)
  cat("  ==========\n")
  dMat <- gw.dist(dp.locat=dp.locat, p=p.max, theta=theta0, longlat=longlat)
  if(!bw.given)
        bw.pmax<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
  res1 <- gwr.aic1(bw.pmax, x, y, kernel,adaptive, dp.locat, dMat)
  AIC.pmax <- res1[[1]]
  if(adaptive)
    cat("    Adaptive bandwidth (number of nearest neighbours): ",bw.pmax, ", AICc value: ", AIC.pmax, "\n")
  else
    cat("    Fixed bandwidth: ",bw.pmax, ", AICc value: ", AIC.pmax, "\n")
  AIC.pi <- 0
  AIC.sel <- 0 
  i <- 0
  p.i <- p.min
  while(p.i<p.max && abs(AIC.pmax-AIC.pi)>=drop.tol)
  {
     cat("  ==========  Try Minkovski distance with p: ", p.i)
     cat("  ==========\n")
     dMat <- gw.dist(dp.locat=dp.locat, p=p.i, theta=theta0, longlat=longlat)
     if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
     bws <- c(bws, bw)
     res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
     AIC.pi <- res1[[1]]
     if(adaptive)
        cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", AIC.pi, "\n")
     else
        cat("    Fixed bandwidth: ",bw, ", AICc value: ", AIC.pi, "\n")
     AICc.vals <- c(AICc.vals, AIC.pi)
     p.vals <- c(p.vals, p.i)
     if(abs(AIC.pi-AIC.sel)>=drop.tol && abs(AIC.pmax-AIC.pi)>=drop.tol)
     {
        AIC.sel <- AIC.pi
        p.dropped <- c(p.dropped, FALSE)
     }
     else
        p.dropped <- c(p.dropped, TRUE)
     p.i <- p.i + interval
  }
  p.vals <- c(p.vals, p.max)
  bws <-c(bws,bw.pmax)
  AICc.vals <- c(AICc.vals, AIC.pmax)
  p.dropped <- c(p.dropped, FALSE)
  res <- list(p.vals=p.vals, cretion.vals=AICc.vals, p.dropped=p.dropped)
  class(res) <- "pvlas"
  res   
}

gwr.mink.pval.backward <- function(formula, data, bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", p.max=2,p.min=0.1,
                       interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given regression data must be Spatial*DataFrame")
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  var.n <- ncol(x)
  longlat <- F
  bw.given <- T
  if(missing(bw))
     bw.given <- F
   
  #############
  # the AICc value calculated from the maximum value of p
  AIC.pmin <- 0
  p.vals <- c(p.min)
  AICc.vals <- c()
  p.dropped <- c()
  bws <- c()
  cat("  ==========  Try MinkoWski distance with the minimum p: ", p.min)
  cat("  ==========\n")
  dMat <- gw.dist(dp.locat=dp.locat, p=p.min, theta=theta0, longlat=longlat)
  if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
  res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
  AIC.pmin <- res1[[1]]
  bws <- c(bw)
  AICc.vals <- c(AIC.pmin)
  p.dropped <- c(FALSE)
  if(adaptive)
    cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", AIC.pmin, "\n")
  else
    cat("    Fixed bandwidth: ",bw, ", AICc value: ", AIC.pmin, "\n")
  AICc.vals <- c(AIC.pmin)
  AIC.pi <- 0
  AIC.sel <- 0 
  i <- 0
  p.i <- p.max
  while(p.i>p.min && abs(AIC.pmin-AIC.pi)>=drop.tol)
  {
     cat("  ==========  Try Minkovski distance with p: ", p.i)
     cat("  ==========\n")
     dMat <- gw.dist(dp.locat=dp.locat, p=p.i, theta=theta0, longlat=longlat)
     if(!bw.given)
        bw<-bw.gwr1(x, y, dp.locat,approach=bw.sel.approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower)
     res1 <- gwr.aic1(bw, x, y, kernel,adaptive, dp.locat, dMat)
     bws<- c(bws, bw)
     AIC.pi <- res1[[1]]
     if(adaptive)
       cat("    Adaptive bandwidth (number of nearest neighbours): ",bw, ", AICc value: ", AIC.pi, "\n")
     else
       cat("    Fixed bandwidth: ",bw, ", AICc value: ", AIC.pi, "\n")
     AICc.vals <- c(AICc.vals, AIC.pi)
     p.vals <- c(p.vals, p.i)
     if(abs(AIC.pi-AIC.sel)>=drop.tol && abs(AIC.pmin-AIC.pi)>=drop.tol)
     {
        AIC.sel <- AIC.pi
        p.dropped <- c(p.dropped, FALSE)
     }
     else
        p.dropped <- c(p.dropped, TRUE)
     p.i <- p.i - interval
  }
  res <- list(p.vals=p.vals, cretion.vals=AICc.vals, p.dropped=p.dropped)
  class(res) <- "pvlas"
  res   
}





















 
