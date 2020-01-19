\name{gwr.collin.diagno}
\alias{gwr.collin.diagno}
\title{Local collinearity diagnostics for basic GWR}
\description{
This function provides a series of local collinearity diagnostics for the independent variables of a basic GWR model.
}
\usage{
gwr.collin.diagno(formula, data, bw, kernel="bisquare",
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, probably calculated 
            by \code{bw.gwr} or \code{bw.gwr.lcr}; fixed (distance) or 
            adaptive bandwidth (number of nearest neighbours)}
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
  \item{corr.mat}{Local correlation matrix}
  \item{VIF}{Local Variance inflation factors (VIFs) matrix}
  \item{local_CN}{Local condition numbers}
  \item{VDP}{Local variance-decomposition proportions}
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with VIF, local_CN, VDP and corr.mat}
}
\references{
Wheeler D, Tiefelsdorf M (2005) Multicollinearity and correlation among local 
regression coefficients in geographically weighted regression. Journal of 
Geographical Systems 7:161-187

Wheeler D (2007) Diagnostic tools and a remedial method for collinearity in 
geographically weighted regression. Environment and Planning A 39:2464-2481

Gollini I, Lu B, Charlton M, Brunsdon C, Harris P (2015) GWmodel: an R Package 
for exploring Spatial Heterogeneity using Geographically Weighted Models. Journal 
of Statistical Software, 63(17):1-50
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR}
\concept{collinearity diagnostics}
