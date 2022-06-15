\name{st.dist}
\alias{st.dist}
\title{Spatio-temporal distance matrix calculation}
\description{
  Calculate a distance vector(matrix) between any GW model calibration point(s) and the data points.
}
\usage{
st.dist(dp.locat, rp.locat, obs.tv, reg.tv,focus=0, p=2, 
        theta=0, longlat=F,lamda=0.05,t.units = "auto",
        ksi=0, s.dMat,t.dMat)
}
\arguments{
  \item{dp.locat}{a numeric matrix of two columns giving the coordinates of the data points}
  \item{rp.locat}{a numeric matrix of two columns giving the coordinates of the GW model calibration points}
  \item{obs.tv}{a vector of time tags for each observation, which could be numeric or of \link{POSIXlt} class}
  \item{reg.tv}{a vector of time tags for each regression location, which could be numeric or of \link{POSIXlt} class}
  \item{focus}{an integer, indexing to the current GW model point, if focus=0, all the distances between all the GW model calibration points and data points will be calculated and a distance matrix will be returned; if 0<focus<length(rp.locat), then the distances between the 'focus'th GW model points and data points will be calculated and a distance vector will be returned}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{lamda}{an parameter between 0 and 1 for calculating spatio-temporal distance} 
  \item{t.units}{	character string to define time unit}
  \item{ksi}{an parameter between 0 and PI for calculating spatio-temporal distance, 
             see details in Wu et al. (2014)}
  \item{s.dMat}{a predifined spatial distance matrix for calculating spatio-temporal distances}
  \item{t.dMat}{a predifined temporal distance matrix for calculating spatio-temporal distances}
}

\value{
 Returns a numeric spatio-temporal distance matrix or vector;
 or a matrix with its rows corresponding to the observations and its columns corresponds to the calibration points.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{gtwr}