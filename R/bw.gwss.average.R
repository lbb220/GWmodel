bw.gwss.average<- function (data, summary.locat, vars, kernel = "bisquare", adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat) 
{
    if (is(data, "Spatial")) 
    {
        p4s <- proj4string(data)
        dp.locat <- coordinates(data)
    }
    else if (is(data, "data.frame") && (!missing(dMat))) 
        data <- data
    else 
        stop("Given data must be a Spatial*DataFrame or data.frame object")
    if (missing(summary.locat)) 
    {
        sp.given <- FALSE
        summary.locat <- data
        sp.locat <- coordinates(summary.locat)
    }
    else 
    {
        sp.given <- T
        if (is(summary.locat, "Spatial")) 
            sp.locat <- coordinates(summary.locat)
        else 
        {
            warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
            summary.locat <- sp.locat
        }
    }
    data <- as(data, "data.frame")
    dp.n <- nrow(data)
    sp.n <- nrow(sp.locat)
    if (missing(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != sp.n) 
            stop("Dimensions of dMat are not correct")
    }
    if (missing(vars)) 
        stop("Variables input error")

    if (missing(dMat)) {
        DM.given <- F
        if (dp.n + dp.n <= 10000) {
            dMat <- gw.dist(dp.locat = dp.locat, rp.locat = dp.locat, 
                p = p, theta = theta, longlat = longlat)
            DM.given <- T
        }
    }
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    if (adaptive) {
        upper <- dp.n
        lower <- 20
    }
    else {
        if (DM.given) {
            upper <- range(dMat)[2]
            lower <- upper/5000
        }
        else {
            dMat <- NULL
            if (p == 2) {
                b.box <- bbox(dp.locat)
                upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 + 
                  (b.box[2, 2] - b.box[2, 1])^2)
                lower <- upper/5000
            }
            else {
                upper <- 0
                for (i in 1:dp.n) {
                  dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                    p = p, theta = theta, longlat = longlat)
                  upper <- max(upper, range(dist.vi)[2])
                }
                lower <- upper/5000
            }
        }
    }
    
   col.nm <- colnames(data)
    var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
    if (length(var.idx) == 0) 
        stop("Variables input doesn't match with data")
    X <- data[, var.idx]
    X <- as.matrix(X)
    colnames(X)<-vars
    var.n <- length(vars)    	
    	bws<-matrix(numeric(2*var.n),nrow=2)
		rownames(bws)<-c('Local Mean bw','Local Median bw')
    colnames(bws) <- colnames(X)
    for(k in 1:var.n)
    {
        bws[1,k] <- gold(gw.mean.cv, lower, upper, adapt.bw = adaptive, 
            X[,k], kernel, adaptive, dp.locat, p, theta, longlat, dMat)
        bws[2,k] <- gold(gw.median.cv, lower, upper, adapt.bw = adaptive, 
            X[,k], kernel, adaptive, dp.locat, p, theta, longlat, dMat)
    } 
    bws
}