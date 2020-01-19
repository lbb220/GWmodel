\name{gwr.mink.approach}
\alias{gwr.mink.approach}
\alias{mink.approach}
\alias{bw.gwr1}
\alias{gwr.aic1}
\alias{gwr.cv1}
\title{Minkovski approach for GWR}
\description{
This function implements the Minkovski approach to select an 'optimum' distance metric
for calibrating a GWR model.
}
\usage{
gwr.mink.approach(formula, data, criterion="AIC", bw, bw.sel.approach = "AIC",adaptive=F, 
              kernel="bisquare", p.vals=seq(from=0.25, to=8, length.out=32), p.inf = T,
                          theta.vals = seq(from=0, to=0.5*pi, length.out=10), verbose=F, 
                          nlower = 10)}

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
  \item{p.vals}{ a collection of positive numbers used as the power of the Minkowski distance}
  \item{p.inf}{if TRUE, Chebyshev distance is tried for model calibration, i.e. p is infinity}
  \item{theta.vals}{a collection of values used as angles in radians to rotate the coordinate system}
  \item{verbose}{if TRUE and bandwidth selection is undertaken, the bandwidth searches are reported}
  \item{nlower}{the minmum number of nearest neighbours if an adaptive kernel is used}
}
\value{
A list of:
  \item{diag.df}{a data frame with four columns (p, theta, bandwidth, AICc/CV), each row corresponds to a calibration}
  \item{coefs.all}{a list class object including all the estimated coefficients}
}
\note{
The function \dQuote{mink.approach} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwr.mink.approach}, while the old name is still kept valid.
}
\references{
Lu, B, Charlton, M, Brunsdon, C & Harris, P(2016). The Minkowski approach for choosing the distance metric in Geographically Weighted Regression. International Journal of Geographical Information Science, 30(2): 351-368.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR}
\concept{Minkowski approach}

