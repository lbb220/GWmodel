\name{gwr.psdm}
\alias{gwr.psdm}
\alias{print.psdmgwr}
\alias{gwr.backfit}
\alias{bw.gwr2}
\title{GWR with parameter-specific distance metrics}
\description{
This function implements GWR with parameter-specific distance metrics to detect variations in
regression relationships across different spatial scales.
}
\usage{
gwr.psdm(formula, data, kernel="bisquare", adaptive=FALSE, criterion="dCVR", 
          max.iterations=1000,threshold=0.000001, dMats,p.vals, theta.vals,
          longlat=FALSE, bws0, bw.seled=rep(F, length(bws0)), approach = "AIC", 
          bws.thresholds=rep(1, length(dMats)), verbose=F, nlower = 10)
\method{print}{psdmgwr}(x, \dots)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{criterion}{criterion for determining the convergence of the back-fitting procedure, could be "CVR" or "dCVR", which corespond to the changing value of RSS (CVR) and the differential version (dCVR), respectively; and "dCVR" is used as default.}
  \item{max.iterations}{maximum number of iterations in the back-fitting procedure}
  \item{threshold}{threshold value to terminate the back-fitting iterations}
  \item{dMats}{a \link{list} of distance matrices used for estimating each specific parameter}
  \item{p.vals}{ a collection of positive numbers used as the power of the Minkowski distance}
  \item{theta.vals}{a collection of values used as angles in radians to rotate the coordinate system}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{bws0}{a \link{vector} of initializing bandwidths for the back-fitting procedure, of which the length should equal to the number of paramters if specified}
  \item{bw.seled}{a \link{vector} of boolean variables to determine whether the corresponding bandwidth should be re-selected or not: if TRUE, the corresponding bandwiths for the specific parameters are                      supposed to be given in bws0; otherwise, the bandwidths for the specific parameters will be selected within the back-fitting iterations.}
  \item{approach}{specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
  \item{bws.thresholds}{threshold values to define whether the bandwidth for a specific parameter has converged or not}
  \item{verbose}{if TRUE and bandwidth selection is undertaken, the bandwidth searches are reported}
  \item{nlower}{the minmum number of nearest neighbours if an adaptive kernel is used}
  \item{x}{an object of class \dQuote{psdmgwr}, returned by the function \link{gwr.psdm}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{psdmgwr}:
 \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with data locations,coefficient estimates from the PSDM GWR model,predicted y values,residuals, coefficient standard errors and t-values in its "data" slot.}
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{bws.vars}{bandwidths used for all the parameters within the back-fitting procedure}
  \item{AICc.vals}{AICc values calculated within the back-fitting procedure}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\note{
This function calibrates a GWR model with parameter-specific distance metrics, 
which by construction also implements a GWR model with parameter-specific bandwidths.
Thus, GWR with flexible bandwidths (FBGWR) and mixed GWR are both special cases 
of PSDM GWR (as are basic GWR and the global regression). Specifically, an FBGWR 
model will be calibrated if no \dQuote{dMats} and \dQuote{p.vals} are specified; a mixed GWR model
will be calibrated if an infinite bandwidth and another regular bandwidth are used for estimating
the global and local parameters (again when no  \dQuote{dMats} and  \dQuote{p.vals} are specified).  
In other words, the \link{gwr.psdm} function is specified with Euclidean distances in both cases.
Note that the results from this function for a mixed GWR model and \link{gwr.mixed} 
might be different, as a back-fitting algorithm is used in \link{gwr.psdm}, while an approximating algorithm is
applied in \link{gwr.mixed}. The gwr.mixed performs better in computational efficiency, 
but poorer in prediction accuracy.
}

\references{

Yang, W. (2014). An Extension of Geographically Weighted Regression with Flexible Bandwidths. 
St Andrews, St Andrews, UK.

Lu, B., Harris, P., Charlton, M., & Brunsdon, C. (2015). Calibrating a Geographically 
Weighted Regression Model with Parameter-specific Distance Metrics. Procedia Environmental Sciences, 26, 109-114.

Lu, B., Brunsdon, C., Charlton, M., & Harris, P. (2017). Geographically weighted 
regression with parameter-specific distance metrics. International Journal of Geographical Information Science, 31, 982-998.

}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
EUDM <- gw.dist(coordinates(londonhp))
#No bandwidth is selected, and bws0 values are used
\dontrun{
res1<- gwr.psdm(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="CVR",kernel="gaussian", 
adaptive=T, bws0=c(100, 100, 100),bw.seled=rep(T, 3), dMats=list(EUDM,EUDM,EUDM))
#FBGWR
res1<- gwr.psdm(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="dCVR",kernel="gaussian")
#Mixed GWR
res3<- gwr.psdm(PURCHASE~FLOORSZ+PROF, data=londonhp, bws0=c(Inf, 100, 100, Inf),
               bw.seled=rep(T, 3),kernel="gaussian", dMats=list(EUDM,EUDM,EUDM))
#PSDM GWR
res4<- gwr.psdm(PURCHASE~FLOORSZ+PROF, data=londonhp, kernel="gaussian", p.vals=c(1,2,3))
}
}
\keyword{multi-scale, flexible bandwidth, parameter-specific distance metrics, GWR}

