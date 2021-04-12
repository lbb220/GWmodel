\name{gwr.lcr}
\alias{gwr.lcr}
\alias{ridge.lm}
\alias{print.gwrlcr}
\title{GWR with a locally-compensated ridge term}
\description{
To address possible local collinearity problems in basic GWR, GWR-LCR finds local ridge parameters at affected
locations (set by a user-specified threshold for the design matrix condition number).
}
\usage{
gwr.lcr(formula, data, regression.points, bw, kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,cv=T,dMat)
\method{print}{gwrlcr}(x, \dots)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame 
                           or SpatialPolygonsDataFrame as defined in package \pkg{sp},
                           or a two-column numeric array}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \code{bw.gwr.lcr};
            fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{lambda}{option for a globally-defined (constant) ridge parameter. Default is lambda=0, which gives a basic GWR fit}
  \item{lambda.adjust}{a locally-varying ridge parameter. Default FALSE, refers to: (i) a basic GWR without
a local ridge adjustment (i.e. lambda=0, everywhere); or (ii) a penalised GWR with a global ridge adjustment 
(i.e. lambda is user-specified as some constant, other than 0 everywhere); if TRUE, use cn.tresh
to set the maximum condition number. Here for locations with a condition number (for its local design matrix) 
above this user-specified threshold, a local ridge parameter is found}
  \item{cn.thresh}{maximum value for condition number, commonly set between 20 and 30}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{cv}{if TRUE, 'cross-validation data will be calculated and returned in the output Spatial*DataFrame}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{x}{an object of class \dQuote{gwrlcr}, returned by the function \link{gwr.lcr}}
  \item{...}{arguments passed through (unused)}
  
}
\value{
A list of class \dQuote{rgwr}:
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object 
            (see package "sp") with coordinates of regression.points in its "data" slot.}
  \item{GW.arguments}{parameters used for the LCR-GWR calibration}
  \item{GW.diagnostic}{diagnostic information is given when data points are also used as regression locations}
  \item{timings}{timing information for running this function}
  \item{this.call}{the function call used.}
}
\references{
Wheeler D (2007) Diagnostic tools and a remedial method for collinearity in geographically weighted regression. Environment and Planning A 39:2464-2481

Brunsdon C, Charlton M, Harris P (2012) Living with collinearity in Local Regression Models. GISRUK 2012, Lancaster, UK

Brunsdon C, Charlton M, Harris P (2012) Living with collinearity in Local Regression Models. Spatial Accuracy 2012, Brazil

Gollini I, Lu B, Charlton M, Brunsdon C, Harris P (2015) GWmodel: an R Package for 
exploring Spatial Heterogeneity using Geographically Weighted Models.  Journal of 
Statistical Software 63(17): 1-50

}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(DubVoter)
require(RColorBrewer)

# Function to find the global condition number (CN)
BKW_cn <- function (X) {
  p <- dim(X)[2]
  Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  Xsvd <- svd(Xscale)$d
  cn <- Xsvd[1] / Xsvd[p]
  cn
}
#
X <- cbind(1,Dub.voter@data[,3:10])
head(X)
CN.global <- BKW_cn(X)
CN.global
\dontrun{
# gwr.lcr function with a global bandwidth to check that the global CN is found
gwr.lcr1 <- gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64, data=Dub.voter, bw=10000000000)
summary(gwr.lcr1$SDF$Local_CN)

# Find and map the local CNs from a basic GWR fit using the lcr-gwr function 
#(note this is NOT the locally-compensated ridge GWR fit as would need to set 
#lambda.adjust=TRUE and cn.thresh=30, say)

bw.lcr2 <- bw.gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64, data=Dub.voter, kernel="bisquare", adaptive=TRUE)
gwr.lcr2 <- gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64, data=Dub.voter, bw=bw.lcr2, kernel="bisquare", adaptive=TRUE)
if(require("RColorBrewer"))
  spplot(gwr.lcr2$SDF,"Local_CN",col.regions=brewer.pal(9,"YlOrRd"),cuts=8,
  main="Local CN")
}
}
\keyword{GWR-LCR}

