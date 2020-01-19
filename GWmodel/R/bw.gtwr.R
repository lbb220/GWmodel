
###Select the bandwidth for GTWR
# Optimize the bandwidth only via the CV or AICc approach
bw.gtwr<-function(formula, data, obs.tv, approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,verbose=T)
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
  
  ####################################Joke
  if(dp.n>3000)
  {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("          -----A kind suggestion from GWmodel development group\n")
   }  
  #################### Recommond to specify a distance matrix
  if (missing(st.dMat))
  {
      DM.given<-F
      if(dp.n + dp.n <= 20000)
      {
        st.dMat <- st.dist(dp.locat, obs.tv=obs.tv, p=p, theta=theta, longlat=longlat,lamda=lamda,t.units = t.units,ksi=ksi)
        DM.given<-T
      }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(st.dMat)
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
      upper<-range(st.dMat[which(!is.infinite(st.dMat))])[2]
      lower<-upper/5000
    }
    else
    {
      st.dMat <- NULL
	    b.box<-bbox(dp.locat)
	    t.intev <- range(obs.tv)
      upper <- max(st.dist(b.box, obs.tv=t.intev, p=p, theta=theta, longlat=longlat,lamda=lamda,t.units = t.units,ksi=ksi))
      lower<-upper/5000
    }
     
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA
    if(approach=="cv"||approach=="CV")
       bw <- gold(gtwr.cv,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat, obs.tv, p, theta, longlat,lamda, t.units, ksi, st.dMat, verbose)
    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
       bw<- gold(gtwr.aic,lower,upper,adapt.bw=adaptive,x,y,kernel,adaptive, dp.locat, obs.tv, p, theta, longlat,lamda, t.units, ksi, st.dMat, verbose)    
    bw

}


####Calculate the CV score with a given bandwidth
##Author: Binbin Lu
gtwr.cv<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, obs.tv, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,verbose=T)
{
   dp.n<-length(dp.locat[,1])
   #########ST Distance matrix is given or not

  if (is.null(st.dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(st.dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################CV
  CV<-numeric(dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-st.dMat[,i]
    else
    {
       dist.vi<-st.dist(dp.locat, obs.tv=obs.tv, focus=i,p=p, theta=theta, longlat=F,lamda=longlat,t.units = t.units,ksi=ksi)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    W.i[i]<-0
    gw.resi<- try(gw_reg(X, Y, W.i, FALSE, i))
    if(!inherits(gw.resi, "try-error"))
    {
      yhat.noi<-X[i,]%*%gw.resi[[1]]
      CV[i]<-Y[i]-yhat.noi

    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  if (!any(is.infinite(CV)))
     CV.score<-t(CV) %*% CV
  else
     {
        CV.score<-Inf
     }
  if(verbose)
  {
    if(adaptive)
      cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
    else
      cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
  }
  CV.score
}

####For select bandwidth, lamda and ksi, the three parameters at the same time
## p3=c(bw, lamda, ksi)
#gtwr.cv.3p<-function(p3, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, obs.tv, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,verbose=T)
#{
#   dp.n<-length(dp.locat[,1])
#   #########ST Distance matrix is given or not
#
#  if (is.null(st.dMat))
#      DM.given<-F
#  else
#  {
#    DM.given<-T
#    dim.dMat<-dim(st.dMat)
#    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
#    stop ("Dimensions of dMat are not correct")
#  }
#  ############################################CV
#  CV<-numeric(dp.n)
#  for (i in 1:dp.n)
#  {
#    if (DM.given)
#         dist.vi<-st.dMat[,i]
#    else
#    {
#       dist.vi<-st.dist(dp.locat, obs.tv=obs.tv, focus=i,p=p, theta=theta, longlat=F,lamda=longlat,t.units = t.units,ksi=ksi)
#    }
#    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
#    W.i[i]<-0
#    gw.resi<- try(gw_reg(X, Y, W.i, FALSE, i))
#    if(!inherits(gw.resi, "try-error"))
#    {
#      yhat.noi<-X[i,]%*%gw.resi[[1]]
#      CV[i]<-Y[i]-yhat.noi
#
#    }
#    else
#    {
#      CV[i]<-Inf
#      break
#    }
#  }
#  if (!any(is.infinite(CV)))
#     CV.score<-t(CV) %*% CV
#  else
#     {
#        CV.score<-Inf
#     }
#  if(verbose)
#  {
#    if(adaptive)
#      cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
#    else
#      cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
#  }
#  CV.score
#}


####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
gtwr.aic<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, obs.tv, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,verbose=T)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   #########Distance matrix is given or not

  if (is.null(st.dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(st.dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  S<-matrix(nrow=dp.n,ncol=dp.n)
  betas <-matrix(nrow=dp.n, ncol=var.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-st.dMat[,i]
    else
    {
       dist.vi <- st.dist(dp.locat, obs.tv=obs.tv, focus=i,p=p, theta=theta, longlat=F,lamda=longlat,t.units = t.units,ksi=ksi)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    res<- try(gw_reg(X,Y,W.i,TRUE,i))
    if(!inherits(res, "try-error"))
    {
      S[i,]<-res[[2]]
      betas[i,] <- res[[1]]
    }
    else
    {
      S[i,]<-Inf
      break
    }
  }

  if (!any(is.infinite(S)))
  {
     AICc<-AICc(Y,X,betas, S)
  }
  else
    AICc<-Inf
  if(verbose)
  {
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
    else
      cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  }
  AICc
}