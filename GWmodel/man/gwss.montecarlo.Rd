\name{gwss.montecarlo}
\alias{gwss.montecarlo}
\alias{montecarlo.gwss}
\title{Monte Carlo (randomisation) test for \link{gwss}}
\description{
This function implements Monte Carlo (randomisation) tests for the GW summary 
statistics found in \link{gwss}.
}
\usage{
gwss.montecarlo(data, vars, kernel = "bisquare", 
                adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, 
                dMat, quantile=FALSE,nsim=99) 
}
\arguments{
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{vars}{a vector of variable names to be summarized}
  \item{bw}{bandwidth used in the weighting function}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calulate the adaptive kernel, and bw correspond to the 
         number of nearest neighbours, default is FALSE.}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{quantile}{if TRUE, median, interquartile range, quantile imbalance will be calculated}
  \item{nsim}{default 99, the number of randomisations}
}
\value{
  \item{test}{probability of the test statistics of the GW summary statistics; 
  if p<0.025 or if p>0.975 then the true local summary statistics can be said to 
  be significantly different (at the 0.95 level) to such a local summary statistics found by chance.}         
}
\note{
  The function \dQuote{montecarlo.gwss} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwss.montecarlo}, while the old name is still kept valid.
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Brunsdon C, Fotheringham AS, Charlton ME (2002) Geographically weighted summary statistics - a framework for localised exploratory data analysis. Computers, Environment and Urban Systems 26:501-524

Harris P, Brunsdon C (2010) Exploring spatial variation and spatial relationships in a freshwater acidification critical load data set for Great Britain using geographically weighted summary statistics. Computers & Geosciences 36:54-70
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
test.lss<-gwss.montecarlo(data=londonhp, vars=c("PURCHASE","FLOORSZ"), bw=5000,
          kernel ="gaussian", dMat=DM,nsim=99)
test.lss
}
}
\keyword{GWSS}
\concept{Monte Carlo}

