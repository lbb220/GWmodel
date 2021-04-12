\name{gwpca}
\alias{gwpca}
\alias{robustSvd}
\alias{rwpca}
\alias{wpca}
\alias{wt.median}
\alias{print.gwpca}
\title{GWPCA}
\description{
This function implements basic or robust GWPCA.
}
\usage{
gwpca(data, elocat, vars, k = 2, robust = FALSE, kernel = "bisquare",
                  adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, cv = T, scores=F,
                  dMat)
\method{print}{gwpca}(x, \dots)
}

\arguments{
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{elocat}{a two-column numeric array or Spatial*DataFrame object for providing evaluation locations, 
                       i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame 
                       as defined in package \pkg{sp}}
  \item{vars}{a vector of variable names to be evaluated}
  \item{k}{the number of retained components; k must be less than the number of variables}
  \item{robust}{if TRUE, robust GWPCA will be applied; otherwise basic GWPCA will be applied}
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
  \item{cv}{If TRUE, cross-validation data will be found that are used to calculate the cross-validation score for the specified bandwidth.}
  \item{scores}{if scores = TRUE, the scores of the supplied data on the principal components will be calculated.}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{x}{an object of class \dQuote{gwpca}, returned by the function \code{\link{gwpca}}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{gwpca}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{pca}{an object of class inheriting from \dQuote{princomp}, see \link{princomp}. }
  \item{loadings}{the localised loadings}
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with local proportions of variance for each 
             principle components, cumulative proportion and winning variable for the 1st principle component in its "data" slot.}
  \item{gwpca.scores}{the localised scores of the supplied data on the principal components }
  \item{var}{The local amount of variance accounted for by each component}
  \item{CV}{Vector of cross-validation data}
  \item{timings}{starting and ending time.}
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Harris P, Brunsdon C, Charlton M (2011)
Geographically weighted principal components analysis.
International Journal of Geographical Information Science 25:1717-1736

Harris P, Brunsdon C, Charlton M, Juggins S, Clarke A (2014) Multivariate spatial 
outlier detection using robust geographically weighted methods.  Mathematical 
 Geosciences 46(1) 1-31

Harris P, Clarke A, Juggins S, Brunsdon C, Charlton M (2014) Geographically 
weighted methods and their use in network re-designs for environmental monitoring. 
Stochastic Environmental Research and Risk Assessment 28: 1869-1887

Harris P, Clarke A, Juggins S, Brunsdon C, Charlton M (2015) Enhancements to a 
geographically weighted principal components analysis in the context of an 
application to an environmental data set.  Geographical Analysis 47: 146-172
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
if(require("mvoutlier") && require("RColorBrewer"))
{
  data(bsstop)
  Data.1 <- bsstop[, 1:14]
  colnames(Data.1)
  Data.1.scaled <- scale(as.matrix(Data.1[5:14]))  # standardised data...
  rownames(Data.1.scaled) <- Data.1[, 1]
  #compute principal components:
  pca <- princomp(Data.1.scaled, cor = FALSE, scores = TRUE)  
  # use covariance matrix to match the following...
  pca$loadings
  data(bss.background)
  backdrop <- function() 
   plot(bss.background, asp = 1, type = "l", xaxt = "n", yaxt = "n", 
   xlab = "", ylab = "", bty = "n", col = "grey")
  pc1 <- pca$scores[, 1]
  backdrop()
  points(Data.1$XCOO[pc1 > 0], Data.1$YCOO[pc1 > 0], pch = 16, col = "blue")
  points(Data.1$XCOO[pc1 < 0], Data.1$YCOO[pc1 < 0], pch = 16, col = "red")
  
  #Geographically Weighted PCA and mapping the local loadings
  # Coordinates of the sites
  Coords1 <- as.matrix(cbind(Data.1$XCOO,Data.1$YCOO)) 
  d1s <- SpatialPointsDataFrame(Coords1,as.data.frame(Data.1.scaled))
  pca.gw <- gwpca(d1s,vars=colnames(d1s@data),bw=1000000,k=10)
  local.loadings <- pca.gw$loadings[, , 1]  
  
  # Mapping the winning variable with the highest absolute loading
  # note first component only - would need to explore all components..
  
  lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))]
  df1p = SpatialPointsDataFrame(Coords1, data.frame(lead = lead.item))
  backdrop()
  colour <- brewer.pal(8, "Dark2")[match(df1p$lead, unique(df1p$lead))]
  plot(df1p, pch = 18, col = colour, add = TRUE)
  legend("topleft", as.character(unique(df1p$lead)), pch = 18, col = 
      brewer.pal(8, "Dark2"))
  backdrop()
  
  #Glyph plots give a view of all the local loadings together
  glyph.plot(local.loadings, Coords1, add = TRUE)
  
  #it is not immediately clear how to interpret the glyphs fully, 
  #so inter-actively identify the full loading information using:
  check.components(local.loadings, Coords1)
  
  # GWPCA with an optimal bandwidth
  bw.choice <- bw.gwpca(d1s,vars=colnames(d1s@data),k=2) 
  pca.gw.auto  <- gwpca(d1s,vars=colnames(d1s@data),bw=bw.choice,k=2)
  # note first component only - would need to explore all components..
  local.loadings <- pca.gw.auto$loadings[, , 1]  
  
  lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))]
  df1p = SpatialPointsDataFrame(Coords1, data.frame(lead = lead.item))
  backdrop()
  colour <- brewer.pal(8, "Dark2")[match(df1p$lead, unique(df1p$lead))]
  plot(df1p, pch = 18, col = colour, add = TRUE)
  legend("topleft", as.character(unique(df1p$lead)), pch = 18, 
  col = brewer.pal(8, "Dark2"))
  
  # GWPCPLOT for investigating the raw multivariate data
  gw.pcplot(d1s, vars=colnames(d1s@data),focus=359, bw = bw.choice) 
}
}
}
\keyword{GWPCA}

