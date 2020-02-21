\name{gwr.model.selection}
\alias{gwr.model.selection}
\alias{model.selection.gwr}
\alias{Generate.formula}
\alias{extract.mat}
\alias{AICc_rss}
\alias{AICc}
\alias{gw.fitted} 
\alias{ehat}
\alias{rss}
\title{Model selection for GWR with a given set of independent variables}
\description{
This function selects one GWR model from many alternatives based on the AICc values.
}
\usage{
gwr.model.selection(DeVar=NULL,InDeVars=NULL, data=list(),bw=NULL,approach="CV",
               adaptive=F,kernel="bisquare",dMat=NULL,p=2, theta=0, longlat=F,
               parallel.method=F,parallel.arg=NULL)
}
\arguments{
  \item{DeVar}{dependent variable}
  \item{InDeVars}{a vector of independent variables for model selection}
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}  
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \code{\link{bw.gwr}}}
  \item{approach}{specified by \bold{CV} (\bold{cv}) for cross validation approach or \bold{AIC} (\bold{aic}) 
                  for selecting bandwidth by AICc values}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}} }
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance}
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{parallel.method}{Specified by `FALSE` for serial approach, 
                         by `"omp"` for multi-thread approach implemented via OpenMP, 
                         by `"cluster"` for multi-process approach implemented via `parallel` package,
                         by `"cuda"` for parallel approach implemented via CUDA}
  \item{parallel.arg}{Set the argument for parallel approach. 
                      If `parallel.method` is `FALSE`, there is no need to set its value. 
                      If `parallel.method` is `"omp"`, its value is used to set how many threads should be created (default by cores of *cores of CPU* - 1).
                      If `parallel.method` is `"cluster"`, its value is used to set how many R session should be created (default by cores of *cores of CPU* - 1).
                      If `parallel.method` is `"cuda"`, its value is used to set how many samples is included in one group during the calibration. This value should not be too big to avoid the overflow of GPU memory. }
}

\value{
A list of:
  \item{model.list}{a list of all the tried GWR models consisted of formulas and variables.}
  \item{GWR.df}{a data frame consited of four columns: bandwidth, AIC, AICc, RSS}
}

\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\note{
The algorithm for selecting GWR models consists of the following four steps:

Step 1. Start by calibrating all the possible bivariate GWR models by sequentially regressing
a single independent variable against the dependent variable; 

Step 2. Find the best performing model which produces the minimum AICc value, and permanently include the corresponding independent
variable in subsequent models; 

Step 3. Sequentially introduce a variable from the remaining
group of independent variables to construct new models with the permanently included independent variables,
and determine the next permanently included variable from the best fitting model that has
the minimum AICc value; 

Step 4. Repeat step 3 until all the independent variables are permanently
included in the model.

In this procedure, the independent variables are iteratively included into the model in a "forward" direction.
Note that there is a clear distinction between the different number of involved variables in a
selection step, which can be called model levels.
}
\references{
Lu, B, Charlton, M, Harris, P, Fotheringham, AS (2014) Geographically weighted regression 
with a non-Euclidean distance metric: a case study using hedonic house price data. 
International Journal of Geographical Information Science 28(4): 660-681
}
\seealso{\code{\link{gwr.model.view}}, \code{\link{gwr.model.sort}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GWR}
\concept{model selection}

