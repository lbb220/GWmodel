\name{gwr.mink.pval}
\alias{gwr.mink.pval}
\alias{gwr.mink.pval.forward}
\alias{gwr.mink.pval.backward}
\alias{plot.pvlas}
\title{Select the values of p for the Minkowski approach for GWR}
\description{
These functions implement heuristics to select the values of p from two intervals: (0, 2] in a 'backward' direction and (2, Inf) in a 'forward' direction.
}
\usage{
gwr.mink.pval(formula, data, criterion="AIC", bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", left.interval=0.25,
                       right.interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
gwr.mink.pval.forward(formula, data, bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", p.max=Inf,p.min=2,
                       interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
gwr.mink.pval.backward(formula, data, bw, bw.sel.approach = "AIC",
                       adaptive=F, kernel="bisquare", p.max=2,p.min=0.1,
                       interval=0.5,drop.tol=3, theta0=0,verbose=F,nlower = 10)
\method{plot}{pvlas}(x, \dots)
                       }

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{criterion}{the criterion used for distance metric selection, AICc ("AICc") or cross-validation ("CV") score; default is "AICc"}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{bw.sel.approach}{approach used to seclect an optimum bandwidth for each calibration if no bandwidth (bw) is given; 
                         specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{left.interval}{ the step-size for searching the left interval (0, 2] in a 'backward' direction}
  \item{right.interval}{the step-size for searching the right interval (2, Inf) in a 'forward' direction}
  \item{p.max}{ the maximum value of p}
  \item{p.min}{ the minimum value of p}
  \item{interval}{the step-size for searching the given interval in a 'backward' or 'forward' direction}
  \item{drop.tol}{ an AICc difference threshold to define whether the values of p to be dropped or not}
  \item{theta0}{a fixed rotation angle in radians}
  \item{verbose}{if TRUE and bandwidth selection is undertaken, the bandwidth searches are reported}
  \item{nlower}{the minmum number of nearest neighbours if an adaptive kernel is used}
  \item{x}{an object of class \dQuote{pvlas}, returned by these functions}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of:
  \item{p.vals}{a vector of tried values of p}
  \item{cretion.vals}{a vector of criterion values (AICc or CV) for tried values of p}
  \item{p.dropped}{a vector of boolean to label whether a value of p to be dropped or not: TRUE means to be dropped and FALSE means to be used for the Minkowski approach}
}
\references{
Lu, B, Charlton, M, Brunsdon, C & Harris, P(2016). The Minkowski approach for choosing the distance metric in Geographically Weighted Regression. International Journal of Geographical Information Science, 30(2): 351-368.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR}
\concept{Minkowski approach}
