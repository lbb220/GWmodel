\name{gwss}
\alias{gwss}
\alias{local.corr}
\alias{print.gwss}
\title{Geographically weighted summary statistics (GWSS)}
\description{
This function calculates basic and robust GWSS.  This includes geographically weighted means, standard deviations and skew. Robust alternatives include geographically weighted medians, inter-quartile ranges and quantile imbalances. This function also calculates basic geographically weighted covariances together with basic and robust
geographically weighted correlations.
}
\usage{
gwss(data, summary.locat,vars,kernel="bisquare",adaptive=FALSE, bw,p=2, 
            theta=0, longlat=F,dMat,quantile=FALSE)
\method{print}{gwss}(x, \dots)
}
\arguments{
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{summary.locat}{a Spatial*DataFrame object for providing summary locations, 
                       i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame 
                       as defined in package \pkg{sp}}
  \item{vars}{a vector of variable names to be summarized}
  \item{bw}{bandwidth used in the weighting function}
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
  \item{quantile}{if TRUE, median, interquartile range, quantile imbalance will be calculated}
  \item{x}{an object of class \dQuote{gwss}, returned by the function \link{gwss}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{lss}:
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp})
             with local means,local standard deviations,local variance,
             local skew,local coefficients of variation, local covariances, 
             local correlations (Pearson's), local correlations (Spearman's),
             local medians, local interquartile ranges, local quantile imbalances and coordinates.
             }
  \item{...}{other information for reporting}          
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Brunsdon C, Fotheringham AS, Charlton ME (2002) Geographically weighted summary statistics - a framework for localised exploratory data analysis. Computers, Environment and Urban Systems 26:501-524

Harris P, Clarke A, Juggins S, Brunsdon C, Charlton M (2014) Geographically 
weighted methods and their use in network re-designs for environmental monitoring. 
Stochastic Environmental Research and Risk Assessment 28: 1869-1887
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
data(EWHP)
data(EWOutline)
head(ewhp)
houses.spdf <- SpatialPointsDataFrame(ewhp[, 1:2], ewhp)
localstats1 <- gwss(houses.spdf, vars = c("PurPrice", "FlrArea"), bw = 50000)
head(data.frame(localstats1$SDF))
localstats1
##A function for mapping data
if(require("RColorBrewer"))
{
   quick.map <- function(spdf,var,legend.title,main.title) 
   {
     x <- spdf@data[,var]
     cut.vals <- pretty(x)
     x.cut <- cut(x,cut.vals)
     cut.levels <- levels(x.cut)
     cut.band <- match(x.cut,cut.levels)
     colors <- brewer.pal(length(cut.levels), "YlOrRd")
     colors <- rev(colors)
     par(mar=c(1,1,1,1))
     plot(ewoutline,col="olivedrab",bg="lightblue1")
     title(main.title)
     plot(spdf,add=TRUE,col=colors[cut.band],pch=16)
     legend("topleft",cut.levels,col=colors,pch=16,bty="n",title=legend.title)
  }
  quick.map(localstats1$SDF, "PurPrice_LM", "1000's Uk Pounds", 
  "Geographically Weighted Mean")
  par(mfrow = c(1, 2))
  quick.map(localstats1$SDF, "PurPrice_LSKe", "Skewness Level", "Local Skewness")
  quick.map(localstats1$SDF, "PurPrice_LSD", "1000's Pounds", "Local Standard Deviation")
  #Exploring Non-Stationarity of Relationships
  quick.map(localstats1$SDF, "Corr_PurPrice.FlrArea", expression(rho), 
  "Geographically Weighted Pearson Correlation")
  #Robust, Quantile Based Local Summary Statistics
  localstats2 <- gwss(houses.spdf, vars = c("PurPrice", "FlrArea"), 
  bw = 50000, quantile = TRUE)
  quick.map(localstats2$SDF, "PurPrice_Median", "1000 UK Pounds", 
  "Geographically Weighted Median House Price")
}
}
}
\keyword{GWSS}

