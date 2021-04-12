\name{bw.gwpca}
\alias{bw.gwpca}
\title{Bandwidth selection for Geographically Weighted Principal Components Analysis (GWPCA)}
\description{
A function for automatic bandwidth selection to calibrate a basic or robust GWPCA via a cross-validation approach only
}
\usage{
bw.gwpca(data,vars,k=2, robust=FALSE,kernel="bisquare",adaptive=FALSE,p=2, 
         theta=0, longlat=F,dMat)
}

\arguments{
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{vars}{a vector of variable names to be evaluated}
  \item{k}{the number of retained components, and it must be less than the number of variables}
  \item{robust}{if TRUE, robust GWPCA will be applied; otherwise basic GWPCA will be applied}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
Returns the adaptive or fixed distance bandwidth
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
\references{
Harris P, Clarke A, Juggins S, Brunsdon C, Charlton M (2015)
Enhancements to a geographically weighted principal components analysis in the context of an application to an environmental data set.
Geographical Analysis 47: 146-172
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWPCA}
\concept{bandwidth selection}

