\name{gwpca.montecarlo.1}
\alias{gwpca.montecarlo.1}
\alias{montecarlo.gwpca.1}
\alias{plot.mcsims}
\title{Monte Carlo (randomisation) test for significance of GWPCA eigenvalue variability
for the first component only - option 1}
\description{
This function implements a Monte Carlo (randomisation) test for a basic or robust 
GW PCA with the bandwidth pre-specified  and constant. The test evaluates whether 
the GW eigenvalues vary significantly across space for the first component only.
}
\usage{
gwpca.montecarlo.1(data, bw, vars, k = 2, nsims=99,robust = FALSE, kernel = "bisquare",
                   adaptive = FALSE,  p = 2, theta = 0, longlat = F, dMat)
\method{plot}{mcsims}(x, sname="SD of local eigenvalues from randomisations", \dots)
}

\arguments{
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwpca};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{vars}{a vector of variable names to be evaluated}
  \item{k}{the number of retained components; k must be less than the number of variables}
  \item{nsims}{the number of simulations for MontCarlo test}
  \item{robust}{if TRUE, robust GWPCA will be applied; otherwise basic GWPCA will be applied}
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
  \item{x}{an object of class \dQuote{mcsims}, returned by the function \link{gwpca.montecarlo.1} or \link{gwpca.montecarlo.2}}
  \item{sname}{the label for the observed value on the plot}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of components:
  \item{actual}{the observed standard deviations (SD) of eigenvalues}
  \item{sims}{a vector of the simulated SDs of eigenvalues}
}
\note{
The function \dQuote{montecarlo.gwpca.1} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwpca.montecarlo.1}, while the old name is still kept valid.
}
\references{
Harris P, Brunsdon C, Charlton M (2011)
Geographically weighted principal components analysis.
International Journal of Geographical Information Science 25:1717-1736
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
data(DubVoter)
DM<-gw.dist(dp.locat=coordinates(Dub.voter))
gmc.res<-gwpca.montecarlo.1(data=Dub.voter, vars=c("DiffAdd", "LARent",
"SC1", "Unempl", "LowEduc"), bw=20,dMat=DM,adaptive=TRUE)
gmc.res
plot(gmc.res)
}
}
\keyword{GWPCA}
\concept{Monte Carlo}

