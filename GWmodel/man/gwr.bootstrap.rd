\name{gwr.bootstrap}
\alias{gwr.bootstrap}
\alias{print.gwrbsm}
\alias{generate.lm.data}
\alias{parametric.bs}
\alias{parametric.bs.local}
\alias{se.bs}
\alias{bias.bs} 
\alias{ci.bs}
\alias{pval.bs}
\alias{gwrtvar}
\alias{gwrt.mlr}
\alias{gwrt.lag}
\alias{gwrt.err}
\alias{gwrt.sma}
\alias{bw.gwr3}
\title{Bootstrap GWR}
\description{
This function implements bootstrap methods to test for coefficient variability 
found from GWR under model assumptions for each of four null hypotheses: MLR, 
ERR, SMA and LAG models.  Global test statistic results are found, as well local 
observation-specific test results that can be mapped. 
}
\usage{
gwr.bootstrap(formula, data, kernel = "bisquare", approach = "AIC",
                 R = 99, k.nearneigh = 4, adaptive = FALSE, p = 2,
                 theta = 0, longlat = FALSE, dMat, verbose = FALSE,
                 parallel.method = FALSE, parallel.arg = NULL)

\method{print}{gwrbsm}(x, \dots)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{approach}{specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
  \item{R}{number of random samples reapted in the bootstrap procedure}
  \item{k.nearneigh}{number of nearest neighbours concerned in calbrating ERR, SMA and LAG models}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
  \item{verbose}{if TRUE and bandwidth selection is undertaken, the bandwidth searches are reported} 
  \item{x}{an object of class \dQuote{gwrbsm}, returned by the function \link{gwr.bootstrap}}
  \item{parallel.method}{ FALSE as default, and the calibration will be conducted traditionally via the serial technique, 
                         "omp": multi-thread technique with the OpenMP API, 
                         "cluster": multi-process technique with the \pkg{parallel} package,
                         "cuda": parallel computing technique with CUDA}
  \item{parallel.arg}{ if parallel.method is not FALSE, then set the argument by following:
                      if parallel.method is "omp", parallel.arg refers to the number of threads used, and its default value is 
                       the number of cores - 1;
                      if parallel.method is "cluster", parallel.arg refers to the number of R sessions used, and its default value is 
                       the number of cores - 1;
                      if parallel.method is "cuda",  parallel.arg refers to the number of calibrations  included in each group, 
                      but note a too large value may cause the overflow of GPU memory. }
  \item{...}{arguments passed through (unused)} 
}
\value{
A list of class \dQuote{gwrbsm}:
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{results}{modified statistics reported from comparisons between GWR and MLR, ERR, SMA and LAG}
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object
        (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y
value,predicted values, coefficient standard errors and bootstrap p-values in its \dQuote{data} slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\note{
This function implements the bootstrap methods introduced in Harris et al. (2017). 
It provides a global test statistic (the modified one given in Harris et al. 2017)
 and a complementary localised version that can be mapped. The bootstrap methods 
 test for coefficient variability found from GWR under model assumptions for each 
 of four null hypotheses: i) multiple linear regression model (MLR); ii) simultaneous
autoregressive error model (ERR); iii) moving average error model (SMA) and iv) simultaneous
autoregressive lag model (LAG).
}

\references{

Harris, P., Brunsdon, C., Lu, B., Nakaya, T., & Charlton, M. (2017). Introducing 
bootstrap methods to investigate coefficient non-stationarity in spatial regression 
models. Spatial Statistics, 21, 241-261.

}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
\dontrun{
#Example with the Georgia educational attainment data
data(Georgia)
data(GeorgiaCounties)
coords <- cbind(Gedu.df$X, Gedu.df$Y)
Gedu.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
#Make a SpatialPolygonDataFrame
require(RColorBrewer)
gSRDF <- SpatialPolygonsDataFrame(polygons(Gedu.counties), over(Gedu.counties, 
                                  Gedu.spdf),match.ID=T)  
mypalette.1 <- brewer.pal(11,"Spectral")
X11(width=9,height=8)                   
spplot(gSRDF, names(gSRDF)[c(5,7:9)], col.regions=mypalette.1,
cuts=10, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Georgia educational attainment predictor data")))
bsm.res <- gwr.bootstrap(PctBach~PctRural+PctEld+PctFB+PctPov, gSRDF, 
                         R=999, longlat=T)
bsm.res
#local bootstrap tests with respect to: MLR, ERR, SMA and LAG models.
mypalette.local.test <- brewer.pal(10,"Spectral")
X11(width=12,height=16)
spplot(bsm.res$SDF, names(bsm.res$SDF)[14:17], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the MLR model 
                       null hypothesis")))

X11(width=12,height=16)
spplot(bsm.res$SDF, names(bsm.res$SDF)[19:22], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the ERR model 
                       null hypothesis")))
X11(width=12,height=16)
spplot(bsm.res$SDF, names(bsm.res$SDF)[24:27], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the SMA model null
                       hypothesis")))

X11(width=12,height=16)
spplot(bsm.res$SDF, names(bsm.res$SDF)[29:32], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the LAG model null
                       hypothesis")))
################################################################################
#Example with Dublin voter data
data(DubVoter)
X11(width=9,height=8)                   
spplot(Dub.voter, names(Dub.voter)[c(5,7,9,10)], col.regions=mypalette.1,
cuts=10, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Dublin voter turnout predictor data")))
bsm.res1 <- gwr.bootstrap(GenEl2004~LARent+Unempl+Age18_24+Age25_44, Dub.voter
                         , R=999)
bsm.res1

#local bootstrap tests with respect to: MLR, ERR, SMA and LAG models.
X11(width=11,height=8)
spplot(bsm.res1$SDF, names(bsm.res1$SDF)[14:17], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the MLR model null
                        hypothesis")))
X11(width=11,height=8)
spplot(bsm.res1$SDF, names(bsm.res1$SDF)[19:22], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the ERR model null
                        hypothesis")))
X11(width=11,height=8)
spplot(bsm.res1$SDF, names(bsm.res1$SDF)[24:27], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the SMA model 
                            null hypothesis")))
X11(width=11,height=8)
spplot(bsm.res1$SDF, names(bsm.res1$SDF)[29:32], col.regions=mypalette.local.test,
cuts=9, par.settings=list(fontsize=list(text=15)),
main=expression(paste("Local p-values for each coefficient of the LAG model 
                            null hypothesis")))
}
}
\keyword{GWR}
\concept{bootstrap method}

