\name{gwr.basic}
\alias{gwr.basic}
\alias{gw_reg}
\alias{gwr_diag}
\alias{Ci_mat}
\alias{F1234.test}
\alias{print.gwrm}
\title{Basic GWR model}
\description{
This function implements basic GWR
}
\usage{
gwr.basic(formula, data, regression.points, bw, kernel="bisquare",
adaptive=FALSE, p=2, theta=0, longlat=F,dMat,F123.test=F,cv=F, W.vect=NULL)
\method{print}{gwrm}(x, \dots)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}; Note that no diagnostic information will returned if it is assigned}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
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
  \item{F123.test}{If TRUE, conduct three seperate F-tests according to Leung et al. (2000).}
  \item{cv}{if TRUE, cross-validation data will be calculated and returned in the output Spatial*DataFrame}
  \item{W.vect}{default NULL, if given it will be used to weight the distance weighting matrix}
  \item{x}{an object of class \dQuote{gwrm}, returned by the function \code{\link{gwr.basic}}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{gwrm}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
  \item{Ftest.res}{results of Leung's F tests when F123.test is TRUE.}
}
\references{
Brunsdon, C, Fotheringham, S, Charlton, M (1996), Geographically Weighted Regression: A Method for Exploring Spatial Nonstationarity. 
Geographical Analysis 28(4):281-298

Charlton, M, Fotheringham, S, and Brunsdon, C (2007), GWR3.0, \url{http://gwr.nuim.ie/}.

Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Leung, Y, Mei, CL, and Zhang, WX (2000), Statistical tests for spatial nonstationarity 
based on the geographically weighted regression model. Environment and Planning A, 32, 9-32.

Lu, B, Charlton, M, Harris, P, Fotheringham, AS (2014) Geographically weighted regression 
with a non-Euclidean distance metric: a case study using hedonic house price data. 
International Journal of Geographical Information Science 28(4): 660-681
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
##Compare the time consumed with and without a specified distance matrix
\dontrun{
system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
            kernel = "gaussian"))
system.time(DM<-gw.dist(dp.locat=coordinates(londonhp)))
system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
            kernel = "gaussian", dMat=DM))

## specify an optimum bandwidth by cross-validation appraoch
bw1<-bw.gwr(PURCHASE~FLOORSZ, data=londonhp, kernel = "gaussian",dMat=DM)
gwr.res1<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw1,kernel = "gaussian", 
                   dMat=DM)
gwr.res1 }
data(LondonBorough)

nsa = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(561900,200900), 
scale = 500, col=1)
\dontrun{
if(require("RColorBrewer"))
{
  mypalette<-brewer.pal(6,"Spectral")
  x11()
  spplot(gwr.res1$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
  ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
  main="GWR estimated coefficients for FLOORSZ with a fixed bandwidth", 
  col.regions=mypalette, sp.layout=list(nsa, londonborough))}
}
\dontrun{
bw2<-bw.gwr(PURCHASE~FLOORSZ,approach="aic",adaptive=TRUE, data=londonhp, 
            kernel = "gaussian", dMat=DM)
gwr.res2<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw2,adaptive=TRUE,
                    kernel = "gaussian", dMat=DM)
gwr.res2
if(require("RColorBrewer"))
{
  x11()
  spplot(gwr.res2$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
  ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
  main="GWR estimated coefficients for FLOORSZ with an adaptive bandwidth", 
  col.regions=mypalette, sp.layout=list(nsa,londonborough))}
}
}
\keyword{GWR}

