\name{gwr.multiscale}
\alias{gwr.multiscale}
\alias{gwr.q2}
\alias{print.multiscalegwr}
\alias{gwr.backfit}
\title{Multiscale GWR}
\description{
This function implements multiscale GWR to detect variations in regression relationships across
different spatial scales. This function can not only find a different bandwidth for each 
relationship but also (and simultaneously) find a different distance metric for 
each relationship (if required to do so).
}
\usage{
gwr.multiscale(formula, data, kernel = "bisquare", adaptive = FALSE,
                 criterion = "dCVR", max.iterations = 2000, threshold =
                 1e-05, dMats, var.dMat.indx, p.vals, theta.vals,
                 longlat = FALSE, bws0, bw.seled, approach = "AIC", bws.thresholds, 
                 bws.reOpts = 5, verbose = F,
                 hatmatrix = T, predictor.centered = rep(T,
                 length(bws0) - 1), nlower = 10, parallel.method = F,
                 parallel.arg = NULL)
\method{print}{multiscalegwr}(x, \dots)
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
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{criterion}{criterion for determining the convergence of the back-fitting procedure, could be "CVR" or "dCVR", which corespond to the changing value of RSS (CVR) and the differential version (dCVR), respectively; and "dCVR" is used as default.}
  \item{max.iterations}{maximum number of iterations in the back-fitting procedure}
  \item{threshold}{threshold value to terminate the back-fitting iterations}
  \item{dMats}{a \link{list} of distance matrices used for estimating each specific parameter}
  \item{var.dMat.indx}{index corresponds to a specific distance matrix for each exploratory variable, if dMats is provided}
  \item{p.vals}{ a collection of positive numbers used as the power of the Minkowski distance}
  \item{theta.vals}{a collection of values used as angles in radians to rotate the coordinate system}
  \item{longlat}{if TRUE, great circle distances will be calculated}
  \item{bws0}{a \link{vector} of initializing bandwidths for the back-fitting procedure, of which the length should equal to the number of paramters if specified}
  \item{bw.seled}{a \link{vector} of boolean variables to determine whether the corresponding bandwidth should be re-selected or not: if TRUE, the corresponding bandwiths for the specific parameters are                      supposed to be given in bws0; otherwise, the bandwidths for the specific parameters will be selected within the back-fitting iterations.}
  \item{approach}{specified by CV for cross-validation approach or by AIC corrected (AICc) approach}
  \item{bws.thresholds}{threshold values to define whether the bandwidth for a specific parameter has converged or not}
  \item{bws.reOpts}{the number times of continually optimizing each parameter-specific bandwidth even though it meets the criterion of convergence, for avoiding  sub-optimal choice due to illusion of                           convergence;}
  \item{verbose}{if TRUE and bandwidth selection is undertaken, the bandwidth searches are reported}
  \item{predictor.centered}{a logical vector of length equalling to the number of predictors, and note intercept is not included; if the element is TRUE, the corresponding predictor will be centered.}
  \item{hatmatrix}{if TRUE the hatmatrix for the whole model will be calculated, and AICc, adjusted-R2 values will be returned accordingly.}
  \item{nlower}{the minmum number of nearest neighbours if an adaptive kernel is used}
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
  \item{x}{an object of class \dQuote{multiscalegwr}, returned by the function \link{gwr.multiscale}}
  \item{...}{arguments passed through (unused)}
}
\value{
A list of class \dQuote{psdmgwr}:
 \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with data locations,coefficient estimates from the PSDM GWR model,predicted y values,residuals, coefficient standard errors and t-values in its "data" slot.}
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{GW.diagnostic}{a list class object including the diagnostic information of the model fitting}
  \item{lm}{an object of class inheriting from \dQuote{lm}, see \link{lm}. }
  \item{bws.vars}{bandwidths used for all the parameters within the back-fitting procedure}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\note{
This function implements multiscale GWR to detect variations in regression 
relationships across different spatial scales. This function can not only find 
a different bandwidth for each relationship, but also (and simultaneously), find
 a different distance metric for each relationship (i.e. Parameter-Specific Distance 
 Metric GWR, i.e. PSDM GWR).  Note that multiscale GWR (MGWR) has also been referred 
 to as flexible bandwidth GWR (FBGWR) and conditional GWR (CGWR) in the literature. 
 All are one and the same model, but where PSDM-GWR additionally provides a different 
 distance metric option for each relationship.  An MGWR model is calibrated if no \dQuote{dMats} 
 and \dQuote{p.vals} are specified; a mixed GWR model will be calibrated if an 
 infinite bandwidth and another regular bandwidth are used for estimating the global and local 
 parameters (again when no \dQuote{dMats} and \dQuote{p.vals} are specified). 
 In other words, the gwr.multiscale function is specified with Euclidean distances 
 in both cases. Note that the results from this function for a mixed GWR model 
 and gwr.mixed might be different, as a back-fitting algorithm is used in
\link{gwr.multiscale}, while an approximating algorithm is applied in gwr.mixed. 
The \link{gwr.mixed} function performs better in computational efficiency, but 
poorer in prediction accuracy.
}

\references{

Yang, W. (2014). An Extension of Geographically Weighted Regression with 
Flexible Bandwidths. St Andrews, St Andrews, UK.

Lu, B., Harris, P., Charlton, M., & Brunsdon, C. (2015). Calibrating a 
Geographically Weighted Regression Model with Parameter-specific Distance 
Metrics. Procedia Environmental Sciences, 26, 109-114.

Lu, B., Brunsdon, C., Charlton, M., & Harris, P. (2017). Geographically weighted 
regression with parameter-specific distance metrics. International Journal of 
Geographical Information Science, 31, 982-998.

Fotheringham, A. S., Yang, W. & Kang, W. (2017). Multiscale Geographically 
Weighted Regression (MGWR). Annals of the American Association of Geographers, 
107, 1247-1265.

Yu, H., A. S. Fotheringham, Z. Li, T. Oshan, W. Kang & L. J. Wolf. 2019. Inference 
in multiscale geographically weighted regression. Geographical Analysis(In press).

Leong, Y.Y., & Yue, J.C. (2017). A modification to geographically weighted 
regression. International Journal of Health Geographics, 16 (1), 11.

Lu, B., Yang, W. Ge, Y. & Harris, P. (2018). Improvements to the calibration of 
a geographically weighted regression with parameter-specific distance metrics 
and bandwidths. Forthcoming Computers, Environment and Urban Systems.

Wolf, L.J, Oshan, T.M, Fotheringham, A.S. (2018). Single and multiscale models of 
process spatial heterogeneity. Geographical Analysis, 50(3): 223-246.

Murakami, D., Lu, B., Harris, P., Brunsdon, C., Charlton, M., Nakaya, T., & Griffith, D. (2019) 
The importance of scale in spatially varying coefficient modelling. 
Forthcoming Annals of the Association of American Geographers.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
EUDM <- gw.dist(coordinates(londonhp))
#No bandwidth is selected, and bws0 values are used
\dontrun{
###Similar as the basic GWR
res1<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="dCVR",kernel="gaussian", 
adaptive=T, bws0=c(100, 100, 100),bw.seled=rep(T, 3), dMats=list(EUDM,EUDM,EUDM))
#FBGWR
res2<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="dCVR",kernel="gaussian",
adaptive=T, bws0=c(100, 100, 100), dMats=list(EUDM,EUDM,EUDM))
#Mixed GWR
res3<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, bws0=c(Inf, 100, 100, Inf),
               bw.seled=rep(T, 3),kernel="gaussian", dMats=list(EUDM,EUDM,EUDM))
#PSDM GWR
res4<- gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, kernel="gaussian", p.vals=c(1,2,3))
}
}
\keyword{multiscale GWR}
\concept{PSDM GWR}
