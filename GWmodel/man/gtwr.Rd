\name{gtwr}
\alias{gtwr}
\alias{print.gtwrm}
\alias{ti.distv}
\alias{ti.distm}
\alias{ti.dist}
\alias{get.before.ti}
\alias{get.ts}
\alias{get.uloat}
\alias{sdist.mat}
\alias{tdist.mat}
\alias{st.dist}
\title{Geographically and Temporally Weighted Regression}
\description{
  A function for calibrating a Geographically and Temporally Weighted Regression 
  (GTWR) model.
}
\usage{
gtwr(formula, data, regression.points, obs.tv, reg.tv, st.bw, kernel="bisquare",
     adaptive=FALSE, p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0,
     st.dMat)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}; Note that no diagnostic information will returned if it is assigned}
  \item{obs.tv}{a vector of time tags for each observation, which could be numeric or of \link{POSIXlt} class}
  \item{reg.tv}{a vector of time tags for each regression location, which could be numeric or of \link{POSIXlt} class}
  \item{st.bw}{spatio-temporal bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
 \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) 
 corresponds to the number of nearest neighbours (i.e. adaptive distance); default 
 is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{lamda}{an parameter between 0 and 1 for calculating spatio-temporal distance} 
  \item{t.units}{	character string to define time unit}
  \item{ksi}{an parameter between 0 and PI for calculating spatio-temporal distance, 
             see details in Wu et al. (2014)} 
  \item{st.dMat}{a pre-specified spatio-temporal distance matrix}
}
\note{
The function implements GTWR model proposed by Huang et al. (2010) and 
Wu et al. (2014).
}

\value{
  A list of class \dQuote{gtwrm}:
  \item{GTW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GTW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,
             GTWR coefficient estimates, y value,predicted values, coefficient standard errors 
             and t-values in its "data" slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\references{
Huang, B., Wu, B., & Barry, M. (2010). Geographically and temporally weighted 
regression for modeling spatio-temporal variation in house prices. International
 Journal of Geographical Information Science, 24, 383-401.
 
Wu, B., Li, R., & Huang, B. (2014). A geographically and temporally weighted 
autoregressive model with application to housing prices. International Journal 
of Geographical Information Science, 28, 1186-1204. 

Fotheringham, A. S., Crespo, R., & Yao, J. (2015). Geographical and Temporal 
Weighted Regression (GTWR). Geographical Analysis, 47, 431-452.

}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}

\keyword{GTWR}

