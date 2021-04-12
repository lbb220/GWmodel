\name{bw.gwr}
\alias{bw.gwr}
\alias{gwr.aic}
\alias{gold}
\alias{AICc_rss1}
\alias{BIC}
\alias{e_vec}
\alias{fitted}
\alias{gw_cv_all}
\alias{gw_cv_all_cuda}
\title{Bandwidth selection for basic GWR}
\description{
  A function for automatic bandwidth selection to calibrate a basic GWR model
}
\usage{
bw.gwr(formula, data, approach="CV", kernel="bisquare",
       adaptive=FALSE, p=2, theta=0, longlat=F, dMat,
       parallel.method=F,parallel.arg=NULL)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{approach}{specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{parallel.method}{ FALSE as default, and the calibration will be conducted traditionally via the serial technique, 
                         "omp": multi-thread technique with the OpenMP API, 
                         "cluster": multi-process technique with the \pkg{parallel} package,
                         "cuda": parallel computing technique with CUDA}
  \item{parallel.arg}{ if parallel.method is not FALSE, then set the argument by following:
                      if parallel.method is "omp", parallel.arg refers to the number of threads used, and its default value is 
                       the number of cores - 1;
                      if parallel.method is "cluster", parallel.arg refers to the number of R sessions used, and its default value is 
                       the number of cores - 1;
                      if parallel.method is "cuda",  parallel.arg refers to the number of calibrations  included in each group, 
                      but note a too large value may cause the overflow of GPU memory. }
}
\note{
For a discontinuous kernel function, a bandwidth can be specified either as a fixed (constant) distance or 
as a fixed (constant) number of local data (i.e. an adaptive distance).  For a continuous kernel function, 
a bandwidth can be specified either as a fixed distance or as a 'fixed quantity that reflects local sample size'  
(i.e. still an 'adaptive' distance but the actual local sample size will be the sample size as functions are continuous).  
In practise a fixed bandwidth suits fairly regular sample configurations whilst an adaptive bandwidth suits highly irregular 
sample configurations. Adaptive bandwidths ensure sufficient (and constant) local information for each local calibration. 
This note is applicable to all GW models
}

\value{
  Returns the adaptive or fixed distance bandwidth
}

\author{Binbin Lu \email{binbinlu@whu.edu.cn}}

\keyword{GWR}
\concept{bandwidth selection}

