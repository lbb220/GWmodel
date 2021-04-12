\name{gwr.lcr.cv.contrib}
\alias{gwr.lcr.cv.contrib}
\title{Cross-validation data at each observation location for the GWR-LCR model}
\description{
This function finds the individual cross-validation score at each observation location, for a GWR-LCR model, for a specified bandwidth.
These data can be mapped to detect unusually high or low cross-validations scores.
}
\usage{
gwr.lcr.cv.contrib(bw,X,Y,locs,kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
}

\arguments{
  \item{bw}{bandwidth used in the weighting function;fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{X}{a numeric matrix of the independent data with an extra column of \dQuote{ones} for the 1st column}
  \item{Y}{a column vector of the dependent data}
  \item{locs}{a two-column numeric array of observation coordinates}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{lambda}{option for a globally-defined (constant) ridge parameter. Default is lambda=0, which gives a basic GWR fit}
  \item{lambda.adjust}{a locally-varying ridge parameter. Default FALSE, refers to: (i) a basic GWR without
a local ridge adjustment (i.e. lambda=0, everywhere); or (ii) a penalised GWR with a global ridge adjustment 
(i.e. lambda is user-specified as some constant, other than 0 everywhere); if TRUE, use cn.tresh
to set the maximum condition number. Here for locations with a condition number (for its local design matrix) 
above this user-specified threshold, a local ridge parameter is found}
  \item{cn.thresh}{maximum value for condition number, commonly set between 20 and 30}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
  \item{CV}{a data vector consisting of squared residuals, whose sum is the cross-validation score for the specified bandwidth.}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR-LCR}
\concept{point-wise cross-validation scores}

