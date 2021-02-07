\name{gw.dist}
\alias{gw.dist}
\alias{gw_dist}
\title{Distance matrix calculation}
\description{
  Calculate a distance vector(matrix) between any GW model calibration point(s) and the data points.
}
\usage{
gw.dist(dp.locat, rp.locat, focus=0, p=2, theta=0, longlat=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp.locat}{a numeric matrix of two columns giving the coordinates of the data points}
  \item{rp.locat}{a numeric matrix of two columns giving the coordinates of the GW model calibration points}
  \item{focus}{an integer, indexing to the current GW model point, if focus=0, all the distances between all the GW model calibration points and data points will be calculated and a distance matrix will be returned; if 0<focus<length(rp.locat), then the distances between the 'focus'th GW model points and data points will be calculated and a distance vector will be returned}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
}

\value{
 Returns a numeric distance matrix or vector;
 matrix with its rows corresponding to the observations and its columns corresponds to the GW model calibration points.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\seealso{\code{\link{dist}} in \pkg{stats}}
\examples{
dp<-cbind(sample(100),sample(100))
rp<-cbind(sample(10),sample(10))
#Euclidean distance metric is used.
dist.v1<-gw.dist(dp.locat=dp, focus=5, p=2, theta=0, longlat=FALSE)
#Manhattan distance metric is used.
#The coordinate system is rotated by an angle 0.5 in radian.
dist.v2<-gw.dist(dp.locat=dp, focus=5, p=1, theta=0.5)
#Great Circle distance metric is used.
dist.v3<-gw.dist(dp.locat=dp, focus=5, longlat=TRUE)
#A generalized Minkowski distance metric is used with p= 0.75 .
#The coordinate system is rotated by an angle 0.8 in radian.
dist.v4<-gw.dist(dp.locat=dp,rp.locat=rp, focus=5, p=0.75,theta=0.8)
################################
#matrix is calculated
#Euclidean distance metric is used.
dist.m1<-gw.dist(dp.locat=dp, p=2, theta=0, longlat=FALSE)
#Manhattan distance metric is used.
#The coordinate system is rotated by an angle 0.5 in radian.
dist.m2<-gw.dist(dp.locat=dp, p=1, theta=0.5)
#Great Circle distance metric is used.
#dist.m3<-gw.dist(dp.locat=dp, longlat=TRUE)
#A generalized Minkowski distance metric is used with p= 0.75 .
#The coordinate system is rotated by an angle 0.8 in radian.
dist.m4<-gw.dist(dp.locat=dp,rp.locat=rp, p=0.75,theta=0.8)
}
\keyword{GW tools}
\concept{distance}
