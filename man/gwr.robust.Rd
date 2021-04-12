\name{gwr.robust}
\alias{gwr.robust}
\title{Robust GWR model}
\description{
This function implements two robust GWR models.
}
\usage{
gwr.robust(formula, data, bw,filtered=FALSE, 
           kernel = "bisquare",adaptive = FALSE, p = 2, 
           theta = 0, longlat = F, dMat, F123.test = F,
           maxiter=20,cut.filter= 3, cut1=2,cut2=3,delta=1.0e-5) }

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{filtered}{default FALSE, the automatic approach is used, if TRUE the filtered data approach is employed, as that described in Fotheringham et al. (2002 p.73-80)}
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
  \item{F123.test}{default FALSE, otherwise calculate F-test results (Leung et al. 2000)}
  \item{maxiter}{default 20, maximum number of iterations for the automatic approach}
  \item{cut.filter}{If filtered is TRUE, it will be used as the residual cutoff for filtering data; default cutoff is 3}
  \item{cut1}{default 2, first cutoff for the residual weighting function. wr(e)=1 if |e|<=cut1*sigma}
  \item{cut2}{default 3, second cutoff for the residual weighting function. wr(e)=(1-(|e|-2)^2)^2 if cut1*sigma<|e|<cut2*sigma, and wr(e)=0 if |e|>=cut2*sigma; cut 1 and cut2 refer to the automatic approach}
  \item{delta}{default 1.0e-5, tolerance of the iterative algorithm}
}
\value{
A list of class \dQuote{gwrm}:
  \item{GW.arguments}{a \link{list} class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a \link{list} class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot. 
             Notably, E_weigts will be also included in the output SDF which represents the residual weighting when automatic approach is used; When the filtered approach is used, E_weight is a vector consisted of 0 and 1, where 0 means outlier to be excluded from calibration.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
  \item{Ftest.res}{results of Leung's F tests when F123.test is TRUE.}
}
\references{

Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Harris P, Fotheringham AS, Juggins S (2010) Robust geographically weighed regression: a technique for quantifying spatial relationships between freshwater acidification critical loads and catchment attributes. Annals of the Association of American Geographers 100(2): 286-306
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
data(DubVoter)
bw.a <- bw.gwr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64,
data=Dub.voter,approach="AICc",kernel="bisquare",adaptive=TRUE)
bw.a
gwr.res <- gwr.basic(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64,
data=Dub.voter,bw=bw.a,kernel="bisquare",adaptive=TRUE,F123.test=TRUE)
print(gwr.res)

# Map of the estimated coefficients for LowEduc
names(gwr.res$SDF)
if(require("RColorBrewer"))
{
  mypalette<-brewer.pal(6,"Spectral")
  X11(width=10,height=12)
  spplot(gwr.res$SDF,"LowEduc",key.space = "right",
  col.regions=mypalette,at=c(-8,-6,-4,-2,0,2,4),
  main="Basic GW regression coefficient estimates for LowEduc")
}
# Robust GW regression and map of the estimated coefficients for LowEduc
rgwr.res <- gwr.robust(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
+Age25_44+Age45_64, data=Dub.voter,bw=bw.a,kernel="bisquare",
adaptive=TRUE,F123.test=TRUE)
print(rgwr.res)
if(require("RColorBrewer"))
{
  X11(width=10,height=12)
  spplot(rgwr.res$SDF, "LowEduc", key.space = "right",
  col.regions=mypalette,at=c(-8,-6,-4,-2,0,2,4),
  main="Robust GW regression coefficient estimates for LowEduc")
}
}
}
\keyword{robust GWR}

