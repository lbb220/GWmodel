gw.average.cv.contrib<- function(bw, X, kernel, adaptive, dp.locat, p, theta, longlat, dMat) 
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
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    	CV.mean <- numeric(dp.n)
	CV.median <- numeric(dp.n)
    for (i in 1:dp.n) {
       	if (DM.given) 
           	dist.vi <- dMat[, i]
        	else {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        	}
		W.i <- matrix(gw.weight(dist.vi, bw, kernel, adaptive),nrow=1)
       	sum.w <- sum(W.i)
       	Wi <- W.i/sum.w
		l.mean<-Wi%*%X
		l.median<- findmedian(X, w = c(Wi))
        Wi <- Wi[-i]
        Wi<-Wi/sum(Wi)
        l.mean.resi <- try(sum(Wi*X[-i]))
        l.median.resi<- try(findmedian(X[-i], w = Wi))
        
        if (!inherits(l.mean.resi, "try-error")) {
            CV.mean[i] <- l.mean - l.mean.resi
        }
        else {
            CV.mean[i] <- Inf
            break
        }
        
        if (!inherits(l.median.resi, "try-error")) {
            CV.median[i] <- l.median - l.median.resi
        }
        else {
            CV.median[i] <- Inf
            break
        }
    }
    CV<-cbind(CV.mean,CV.median)
    colnames(CV)<-c('Local Mean','Local Median')
    CV
}