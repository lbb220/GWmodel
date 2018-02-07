#Calculate the cv score for each column, not the whole X of all the variables
gw.mean.cv <- function (bw, x, kernel, adaptive, dp.locat, p, theta, longlat, 
    dMat) 
{	
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
	  CV.mean <- numeric(dp.n)
    for (i in 1:dp.n) 
    {
       	if (DM.given) 
           	dist.vi <- dMat[, i]
       	else 
        {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        }
       	W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
       	sum.w <- sum(W.i)
       	Wi <- W.i/sum.w
       	l.mean<-sum(Wi*x)
        W.i[i] <- 0 
        W.i<-W.i/sum(W.i)
        l.mean.resi <- try(sum(W.i*x))
        if (!inherits(l.mean.resi, "try-error")) {
            CV.mean[i] <- l.mean - l.mean.resi
        }
        else {
            CV.mean[i] <- Inf
            break
        }
        
    	}
    if (!any(is.infinite(CV.mean))) 
        CV.mean.score <- t(CV.mean) %*% CV.mean
    else {
        CV.mean.score <- Inf
    }

    CV.mean.score
}

gw.median.cv <- function (bw, X, kernel, adaptive, dp.locat, p, theta, longlat, dMat) 
{
	findmedian <- function(x, w) {
        lw <- length(w)
        xo <- sort(x)
        wo <- w[order(x)]
       
        cond <- max({ cumsum(wo) <= 0.5} * seq(1:lw))
            if (cond == 0) 
                cond <- 1
        xo[cond]
    }
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else 
    {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
	CV.median <- numeric(dp.n)
    for (i in 1:dp.n) {
       	if (DM.given) 
           	dist.vi <- dMat[, i]
        	else {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        	}
       	W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
       	sum.w <- sum(W.i)
       	Wi <- W.i/sum.w
		l.median<- findmedian(X, w = Wi)
		Wi <- Wi[-i]
        Wi<-Wi/sum(Wi)
        l.median.resi<- try(findmedian(X[-i], w = Wi))
        if (!inherits(l.median.resi, "try-error")) {
            CV.median[i] <- l.median - l.median.resi
        }
        else {
            CV.median[i] <- Inf
            break
        }
        
    	}
        if (!any(is.infinite(CV.median))) 
        CV.median.score <- t(CV.median) %*% CV.median
    else {
        CV.median.score <- Inf
    }
  
   CV.median.score
  
}

# bws should be a matrix of 2*var.n, variable specific
gw.average.cv <- function (bws, X, kernel, adaptive, dp.locat, p, theta, longlat, 
    dMat) 
{
	findmedian <- function(x, w) {
        lw <- length(w)
        xo <- sort(x)
        wo <- w[order(x)]
       
        cond <- max({ cumsum(wo) <= 0.5} * seq(1:lw))
            if (cond == 0) 
                cond <- 1
        xo[cond]
    }

	
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else 
    {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }

    vars<-colnames(X)
    var.n <- length(vars)
    if(dim(bws)[1]!=2||dim(bws)[2]!=var.n)
      stop("Dimensions of bws are not correct (see the mannual for tips)!") 
	  CV.mean <- matrix(numeric(dp.n*var.n),ncol=var.n)
	  CV.median <- matrix(numeric(dp.n*var.n),ncol=var.n)
    CV.mean.score <- numeric(var.n)
    CV.median.score <- numeric(var.n)
    for(k in 1:var.n)
    {
      x <- X[,k] 
      for (i in 1:dp.n) 
      {
         	if (DM.given) 
             	dist.vi <- dMat[, i]
          	else {
              dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                  p = p, theta = theta, longlat = longlat)
          	}
         	W.i <- gw.weight(dist.vi, bws[1,k], kernel, adaptive)
         	sum.w <- sum(W.i)
         	Wi <- W.i/sum.w
    	    l.mean<-sum(Wi*x)
          W.i[i] <- 0
          W.i<-W.i/sum(W.i)
          l.mean.resi <- try(sum(W.i*x))   
          if (!inherits(l.mean.resi, "try-error")) {
              CV.mean[i,k] <- l.mean - l.mean.resi
          }
          else 
          {
              CV.mean[i,k] <- Inf
              break
          }   
          W.i <- gw.weight(dist.vi, bws[2,k], kernel, adaptive)
         	sum.w <- sum(W.i)   
          Wi <- W.i/sum.w   
    	    l.median<- findmedian(x, w = Wi)
          Wi <- Wi[-i]
          Wi <- Wi/sum(Wi)
          l.median.resi<- try(findmedian(X[-i], w = Wi))
          if (!inherits(l.median.resi, "try-error")) 
          {
              CV.median[i,k] <- l.median - l.median.resi
          }
          else {
              CV.median[i,k] <- Inf
              break
          } 
      	}
       if (!any(is.infinite(CV.mean[,k]))) 
           CV.mean.score[k] <- diag(t(CV.mean[,k]) %*% CV.mean[,k])
       else 
           CV.mean.score[k] <- Inf
       
       if (!any(is.infinite(CV.median[,k]))) 
          CV.median.score[k] <- t(CV.median[,k]) %*% CV.median[,k]
       else 
          CV.median.score[k] <- Inf
    }
  	CV.score<-list(CV.mean.score=CV.mean.score,CV.mean=CV.mean,CV.median.score=CV.median.score,CV.median=CV.median)
   
    CV.score
}