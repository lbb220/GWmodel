\name{bw.gtwr}
\alias{bw.gtwr}
\alias{gtwr.cv}
\alias{gtwr.aic}
\title{Bandwidth selection for GTWR}
\description{
  A function for automatic bandwidth selection to calibrate a GTWR model
}
\usage{
bw.gtwr(formula, data, obs.tv, approach="CV",kernel="bisquare",adaptive=FALSE, 
        p=2, theta=0, longlat=F,lamda=0.05,t.units = "auto",ksi=0, st.dMat,
        verbose=T)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{obs.tv}{a vector of time tags for each observation, which could be numeric or of \link{POSIXlt} class}
  \item{approach}{specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
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
  \item{verbose}{logical variable to define whether show the selection procedure}
}
\note{
The function is developed according to the articles by Huang et al. (2010) and 
Wu et al. (2014).
}

\value{
  Returns the adaptive or fixed distance bandwidth
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
\concept{bandwidth selection}
