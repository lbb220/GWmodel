\name{ggwr.basic}
\alias{ggwr.basic}
\alias{gwr.generalised}
\alias{gwr.binomial}
\alias{gwr.binomial.wt}
\alias{gwr.poisson}
\alias{gwr.poisson.wt}
\alias{gwr.fitted}
\alias{print.ggwrm}
\title{Generalised GWR models with Poisson and Binomial options}
\description{
This function implements generalised GWR
}
\usage{
ggwr.basic(formula, data, regression.points, bw, family =
                 "poisson", kernel = "bisquare", adaptive = FALSE, cv =
                 T, tol = 1e-05, maxiter = 20, p = 2, theta = 0,
                 longlat = F, dMat, dMat1)

 \method{print}{ggwrm}(x, \dots)
 }

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by bw.ggwr();fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{family}{a description of the error distribution and link function to
          be used in the model, which can be specified by \dQuote{poisson} or \dQuote{binomial}}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube:  wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar:   wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{cv}{if TRUE, cross-validation data will be calculated}
  \item{tol}{the threshold that determines the convergence of the IRLS procedure}
  \item{maxiter}{the maximum number of times to try the IRLS procedure}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix between regression points and observations, it can be calculated by the function \code{\link{gw.dist}}}
  \item{dMat1}{a square distance matrix between each pair of observations, it can be calculated by the function \code{\link{gw.dist}}}
  \item{x}{an object of class \dQuote{ggwrm}, returned by the function \link{gwr.generalised}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{ggwrm}:
  \item{GW.arguments}{a \link{list} class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a \link{list} class object including the diagnostic information of the model fitting}
  \item{glm.res}{an object of class inheriting from \dQuote{glm} which inherits from the class \dQuote{lm}, see \link{glm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot.}
  \item{CV}{a data vector consisting of the cross-validation data}
}
\note{
Note that this function calibrates a Generalised GWR model via an approximating 
algorithm, which is different from the back-fitting algorithm used in the GWR4 software
by  Tomoki Nakaya.
}
\references{
Nakaya, T., A. S. Fotheringham, C. Brunsdon & M. Charlton (2005) Geographically 
weighted Poisson regression for disease association mapping. Statistics in Medicine, 
24, 2695-2717.

Nakaya, T., M. Charlton, S. Fotheringham & C. Brunsdon. 2009. How to use SGWRWIN (GWR4.0). 
Maynooth, Ireland: National Centre for Geocomputation.

Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
\dontrun{
DM<-gw.dist(dp.locat=coordinates(londonhp))
bw.f1 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM)
res.poisson<-ggwr.basic(BATH2~FLOORSZ, bw=bw.f1,data=londonhp, dMat=DM)
bw.f2 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM,family ="binomial")
res.binomial<-ggwr.basic(BATH2~FLOORSZ, bw=bw.f2,data=londonhp, dMat=DM,
              family ="binomial")
}
}
\keyword{generalised GWR}
