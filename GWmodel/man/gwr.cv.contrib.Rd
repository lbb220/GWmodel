\name{gwr.cv.contrib}
\alias{gwr.cv.contrib}
\title{Cross-validation data at each observation location for a basic GWR model}
\description{
This function finds the individual cross-validation score at each observation location, for a basic GWR model, for a specified bandwidth.
These data can be mapped to detect unusually high or low cross-validations scores.
}
\usage{
gwr.cv.contrib(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, 
               theta=0, longlat=F,dMat,
               parallel.method=F, parallel.arg=NULL)
}

\arguments{
  \item{bw}{bandwidth used in the weighting function;fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{X}{a numeric matrix of the independent data with an extra column of \dQuote{ones} for the 1st column}
  \item{Y}{a column vector of the dependent data}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{dp.locat}{a two-column numeric array of observation coordinates}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{parallel.method}{Specified by `FALSE` for serial approach, 
                         by `"omp"` for multi-thread approach implemented via OpenMP, 
                         by `"cluster"` for multi-process approach implemented via `parallel` package,
                         by `"cuda"` for parallel approach implemented via CUDA}
  \item{parallel.arg}{Set the argument for parallel approach. 
                      If `parallel.method` is `FALSE`, there is no need to set its value. 
                      If `parallel.method` is `"omp"`, its value is used to set how many threads should be created (default by cores of *cores of CPU* - 1).
                      If `parallel.method` is `"cluster"`, its value is used to set how many R session should be created (default by cores of *cores of CPU* - 1).
                      If `parallel.method` is `"cuda"`, its value is used to set how many samples is included in one group during the calibration. This value should not be too big to avoid the overflow of GPU memory. }
}
\value{
  \item{CV}{a data vector consisting of squared residuals, whose sum is the cross-validation score for the specified bandwidth.}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR}
\concept{point-wise cross-validation scores}
