\name{gwr.mixed}
\alias{gwr.mixed}
\alias{gwr.mixed.2}
\alias{gwr.mixed.trace}
\alias{print.mgwr}
\alias{gwr.q}
\title{Mixed GWR}
\description{
This function implements mixed (semiparametric) GWR
}
\usage{
gwr.mixed(formula, data, regression.points, fixed.vars,
                     intercept.fixed=FALSE, bw, diagnostic=T, kernel="bisquare", 
                     adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{fixed.vars}{independent variables that appeared in the formula that are to be treated as global}
  \item{intercept.fixed}{logical, if TRUE the intercept will be treated as global}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{diagnostic}{logical, if TRUE the diagnostics will be calculated}
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
}
\value{
A list of class \dQuote{mgwr}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{aic}{AICc value from this calibration}
  \item{df.used}{ effective degree of freedom}
  \item{rss}{residual sum of squares}
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with coefficient estimates in its "data" slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Brunsdon C, Fotheringham AS, Charlton ME (1999) Some notes on parametric signficance 
tests for geographically weighted regression. Journal of Regional Science 39(3):497-524

Mei L-M, He S-Y, Fang K-T (2004) A note on the mixed geographically weighted regression 
model. Journal of regional science 44(1):143-157

Mei L-M, Wang N, Zhang W-X (2006) Testing the importance of the explanatory variables 
in a mixed geographically weighted regression model. Environment and Planning A 38:587-598

Nakaya T, Fotheringham AS, Brunsdon C, Charlton M (2005) Geographically Weighted Poisson Regression for Disease Association Mapping,
Statistics in Medicine 24: 2695-2717

Nakaya T et al. (2011) GWR4.0, \url{http://gwr.nuim.ie/}.

}
\note{
For an alternative formulation of mixed GWR, please refer to GWR 4, which provides useful tools for automatic bandwidth selection.
This windows-based software also implements generalised mixed GWR.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{multiscale GWR}
\concept{mixed GWR}

