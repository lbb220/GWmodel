\name{gwda}
\alias{gwda}
\alias{print.gwda}
\alias{grouping.xy}
\alias{wqda}
\alias{wlda}
\alias{splitx}
\alias{wmean}
\alias{wvarcov}
\alias{wprior}
\alias{confusion.matrix}
\title{GW Discriminant Analysis}
\description{
This function implements GW discriminant analysis, where location-wise probabilities and their associated entropy are also calculated.
}
\usage{
gwda(formula, data, predict.data,validation = T, COV.gw=T, 
                 mean.gw=T, prior.gw=T, prior=NULL, wqda =F,
                kernel = "bisquare", adaptive = FALSE, bw,
                 p = 2, theta = 0, longlat = F,dMat)
\method{print}{gwda}(x, \dots)
}

\arguments{
  \item{formula}{Model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame for training, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{predict.data}{a Spatial*DataFrame object for prediction, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as 
  defined in package \pkg{sp}; if it is not given, the traing data will be predicted using leave-one-out cross-validation.}
  \item{validation }{If TRUE, the results from the prediction will be validated and the correct proportion will be calculated.}
  \item{COV.gw}{if true, localised variance-covariance matrix is used for GW discriminant analysis; otherwise, global variance-covariance matrix is used}
  \item{mean.gw}{if true, localised mean is used for GW discriminant analysis; otherwise, global mean is used}
  \item{prior.gw}{if true, localised prior probability is used for GW discriminant analysis; otherwise, fixed prior probability is used}
  \item{prior}{a vector of given prior probability}
  \item{wqda}{if TRUE, weighted quadratic discriminant analysis will be applied; otherwise weighted linear discriminant analysis will be applied}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwpca};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{x}{an object of class \dQuote{gwda}}
  \item{...}{arguments passed through (unused)}
}
\value{
An object of class \dQuote{gwda}. This includes a SpatialPointsDataFrame (may be gridded) or 
SpatialPolygonsDataFrame object, SDF, (see package \dQuote{sp}) with, following the use of new version of \link{gwda}, the probabilities for
each level, the highest probabiliity and the entropy of the probabilities in its \dQuote{data} slot.
}

\references{
Brunsdon, C, Fotheringham S,  and Charlton, M (2007),
Geographically Weighted Discriminant Analysis, Geographical Analysis 39:376-396

Lu B, Harris P, Charlton M, Brunsdon C (2014) The GWmodel R Package: further 
topics for exploring Spatial Heterogeneity using Geographically Weighted Models.
Geo-spatial Information Science 17(2): 85-101
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
 data(USelect)
 dMat <- gw.dist(coordinates(USelect2004))
 bw <- bw.gwda(winner~unemploy+pctcoled+PEROVER65+pcturban+WHITE,data=USelect2004,
 adaptive=TRUE,dMat=dMat)
 ge.gwda <- gwda(winner~unemploy+pctcoled+PEROVER65+pcturban+WHITE,data=USelect2004,
 bw=bw,adaptive=TRUE,dMat=dMat)
 table(USelect2004$winner,ge.gwda$SDF$group.predicted)
 spplot(ge.gwda$SDF, "entropy")
 }
}
\keyword{GWDA}

