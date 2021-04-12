\name{gwr.basic}
\alias{gwr.basic}
\alias{gw_reg}
\alias{gwr_diag}
\alias{Ci_mat}
\alias{gw_local_r2}
\alias{gw_reg_1}
\alias{gw_reg_2}
\alias{gw_reg_all}
\alias{gw_reg_all_cuda}
\alias{gw_cv_all_omp}
\alias{gw_reg_all_omp}
\alias{trhat2}
\alias{F1234.test}
\alias{print.gwrm}
\title{Basic GWR model}
\description{
This function implements basic GWR
}
\usage{
gwr.basic(formula, data, regression.points, bw, kernel="bisquare",
adaptive=FALSE, p=2, theta=0, longlat=F,dMat,F123.test=F,cv=F, W.vect=NULL,
parallel.method=FALSE,parallel.arg=NULL)
\method{print}{gwrm}(x, \dots)
}

\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{regression.points}{a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}; Note that no diagnostic information will returned if it is assigned}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
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
  \item{F123.test}{If TRUE, conduct three seperate F-tests according to Leung et al. (2000).}
  \item{cv}{if TRUE, cross-validation data will be calculated and returned in the output Spatial*DataFrame}
  \item{W.vect}{default NULL, if given it will be used to weight the distance weighting matrix}
  \item{x}{an object of class \dQuote{gwrm}, returned by the function \code{\link{gwr.basic}}}
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
A list of class \dQuote{gwrm}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with fit.points,GWR coefficient estimates, y value,predicted values, coefficient standard errors and t-values in its "data" slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
  \item{Ftest.res}{results of Leung's F tests when F123.test is TRUE.}
}
\references{
Brunsdon, C, Fotheringham, S, Charlton, M (1996), Geographically Weighted Regression: A Method for Exploring Spatial Nonstationarity. 
Geographical Analysis 28(4):281-298

Charlton, M, Fotheringham, S, and Brunsdon, C (2007), GWR3.0, \url{http://gwr.nuim.ie/}.

Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Leung, Y, Mei, CL, and Zhang, WX (2000), Statistical tests for spatial nonstationarity 
based on the geographically weighted regression model. Environment and Planning A, 32, 9-32.

Lu, B, Charlton, M, Harris, P, Fotheringham, AS (2014) Geographically weighted regression 
with a non-Euclidean distance metric: a case study using hedonic house price data. 
International Journal of Geographical Information Science 28(4): 660-681

OpenMP: \url{https://www.openmp.org/}

CUDA: \url{https://developer.nvidia.com/cuda-zone}

R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical
Computing, Vienna, Austria.  \url{https://www.R-project.org/}.
}
\note{
Requirements of using CUDA for high-performence computation in GWR functions:

To run GWR-CUDA (i.e. parallel.method is pecified as \dQuote{cuda}) with  
\code{gwr.basic} , \code{bw.gwr} and \code{gwr.model.selection}, 
the following conditions are required:

1. There is at least one NVIDIA GPU supporting CUDA equipped on user's computer.

2. CUDA (>10.2) are installed with the environment variable `CUDA_HOME` set properly. 

3. The package should re-built from source. 
   - For Linux user, `GWmodelCUDA` will be automatically built if CUDA toolkit could be detected 
     by the complier. 
   - For Windows user, `GWmodelCUDA.dll` and `GWmodelCUDA.lib` will be automatically downloaded;
     however, we would recommend users to build the `GWmodelCUDA` library manually to avoid some potentially
     unknown issues, and an `CMakeLists.txt` file is provided for this procedure.

If any condition above is not satisfied, the GWR-CUDA will not work even though the \dQuote{parallel.method} is 
specified as \dQuote{cuda}.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
##Compare the time consumed with and without a specified distance matrix
\dontrun{
system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
            kernel = "gaussian"))
system.time(DM<-gw.dist(dp.locat=coordinates(londonhp)))
system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
            kernel = "gaussian", dMat=DM))

## specify an optimum bandwidth by cross-validation appraoch
bw1<-bw.gwr(PURCHASE~FLOORSZ, data=londonhp, kernel = "gaussian",dMat=DM)
gwr.res1<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw1,kernel = "gaussian", 
                   dMat=DM)
gwr.res1 }
data(LondonBorough)

nsa = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(561900,200900), 
scale = 500, col=1)
\dontrun{
if(require("RColorBrewer"))
{
  mypalette<-brewer.pal(6,"Spectral")
  x11()
  spplot(gwr.res1$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
  ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
  main="GWR estimated coefficients for FLOORSZ with a fixed bandwidth", 
  col.regions=mypalette, sp.layout=list(nsa, londonborough))}
}
\dontrun{
bw2<-bw.gwr(PURCHASE~FLOORSZ,approach="aic",adaptive=TRUE, data=londonhp, 
            kernel = "gaussian", dMat=DM)
gwr.res2<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw2,adaptive=TRUE,
                    kernel = "gaussian", dMat=DM)
gwr.res2
if(require("RColorBrewer"))
{
  x11()
  spplot(gwr.res2$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
  ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
  main="GWR estimated coefficients for FLOORSZ with an adaptive bandwidth", 
  col.regions=mypalette, sp.layout=list(nsa,londonborough))}
}
\dontrun{
  ############HP-GWR test code
  simulate.data.generator <- function(data.length) {
  x1 <- rnorm(data.length)
  x2 <- rnorm(data.length)
  x3 <- rnorm(data.length)
  lon <- rnorm(data.length, mean = 533200, sd = 10000)
  lat <- rnorm(data.length, mean = 159400, sd = 10000)
  y <- x1 + 5 * x2 + 2.5 * x3 + rnorm(data.length)
  simulate.data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, lon = lon, lat = lat)
  coordinates(simulate.data) <- ~ lon + lat
  names(simulate.data)
  return(simulate.data)
}
simulate.data <- simulate.data.generator(10000)
adaptive = TRUE

## GWR (not parallelized)
bw.CV.s <- bw.gwr(data = simulate.data, formula = y ~ x1 + x2 + x3, approach="CV", 
                  kernel = "gaussian", adaptive = adaptive, parallel.method = FALSE)
model.s <- gwr.model.selection(DeVar = "y", InDeVars = c("x1", "x2", "x3"), data = simulate.data, 
                              bw = bw.CV.s, approach="AIC", kernel = "gaussian", adaptive = T, 
                              parallel.method = FALSE)
system.time(
  betas.s <- gwr.basic(data = simulate.data, formula = y ~ x1 + x2 + x3, bw = bw.CV.s, 
                       kernel = "gaussian", adaptive = TRUE)
)

## GWR-Omp
bw.CV.omp <- bw.gwr(data = simulate.data, formula = y ~ x1 + x2 + x3, approach="CV", 
                    kernel = "gaussian", adaptive = adaptive, parallel.method = "omp")
model.omp <- gwr.model.selection(DeVar = "y", InDeVars = c("x1", "x2", "x3"), data = simulate.data, 
                                bw = bw.CV.omp, approach="AIC", kernel = "gaussian", adaptive = T, 
                                parallel.method = "omp")
system.time(
  betas.omp <- gwr.basic(data = simulate.data, formula = y ~ x1 + x2 + x3, bw = bw.CV.omp, 
                        kernel = "gaussian", adaptive = T, parallel.method = "omp"))

## GWR-CUDA
bw.CV.cuda <- bw.gwr(data = simulate.data, formula = y ~ x1 + x2 + x3, approach="CV", 
                     kernel = "gaussian", adaptive = adaptive, parallel.method = "cuda", 
                     parallel.arg = 6*16)
model.cuda <- gwr.model.selection(DeVar = "y", InDeVars = c("x1", "x2", "x3"), data = simulate.data, 
                                 bw = bw.CV.cuda, approach="AIC", kernel = "gaussian", adaptive = T, 
                                 parallel.method = "cuda", parallel.arg = 6*16)
system.time(
  betas.cuda <- gwr.basic(data = simulate.data, formula = y ~ x1 + x2 + x3, bw = bw.CV.cuda, 
                          kernel = "gaussian", adaptive = T, parallel.method = "cuda", 
                          parallel.arg = 6*8))
}
}
\keyword{GWR}