\name{gwr.hetero}
\alias{gwr.hetero}
\title{Heteroskedastic GWR }
\description{
This function implements a heteroskedastic GWR model
}
\usage{
gwr.hetero(formula, data, regression.points, bw, kernel="bisquare",
                    adaptive=FALSE, tol=0.0001,maxiter=50,verbose=T,
                    p=2, theta=0, longlat=F,dMat)}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{tol}{the threshold that determines the convergence of the iterative procedure}
  \item{maxiter}{the maximum number of times to try the iterative procedure}
  \item{verbose}{logical, if TRUE verbose output will be made from the iterative procedure}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with coefficient estimates in its "data" slot.}
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Harris P, Fotheringham AS, Juggins S (2010) Robust geographically weighed regression: 
a technique for quantifying spatial relationships between freshwater acidification 
critical loads and catchment attributes. Annals of the Association of American Geographers 100(2): 286-306

Harris P, Brunsdon C, Fotheringham AS (2011) Links, comparisons and extensions of the geographically 
weighted regression model when used as a spatial predictor.  Stochastic Environmental Research 
and Risk Assessment 25:123-138

}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{Heteroskedastic GWR}

