\name{gwr.scalable}
\alias{gwr.scalable}
\alias{scgwr_pre}
\alias{scgwr_loocv}
\alias{scgwr_reg}
\alias{AICc1}
\alias{gwr.scalable.loocv}
\alias{gwr_diag1}
\alias{print.scgwrm}
\title{Scalable GWR}
\description{
This function implements Scalable GWR for large dataset
}
\usage{
gwr.scalable(formula, data, bw.adapt=100, kernel = "gaussian", polynomial = 4, 
             p = 2, theta = 0, longlat = F, dMat)
\method{print}{scgwrm}(x, \dots)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw.adapt}{adaptive bandwidth (i.e. number of nearest neighbours) used for geographically weighting}
  \item{kernel}{Kernel function to calculate the spatial weights, but note only two continuous functions available:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);}
  \item{polynomial}{Degree of the polyunomial to approximate the kernel function, and default is 4.}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{x}{an object of class \dQuote{scgwrm}, returned by the function \code{\link{gwr.scalable}}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{scgwrm}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot.}
  \item{timings}{starting and ending time.}
}
\references{
Murakami, D., N. Tsutsumida, T. Yoshida, T. Nakaya & B. Lu (2019) Scalable GWR: A 
linear-time algorithm for large-scale geographically weighted regression with polynomial kernels. arXiv:1905.00266.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
require(spData)
data(boston)
boston <- boston.c
coordinates(boston) <- ~ LON + LAT
res <- gwr.scalable(formula = MEDV ~ CRIM + ZN + INDUS + CHAS + AGE, data = boston, bw.adapt = 100)
res
}
}
\keyword{Scalable GWR}

