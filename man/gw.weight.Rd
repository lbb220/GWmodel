\name{gw.weight}
\alias{gw.weight}
\alias{gw.weight.bis}
\alias{gw.weight.bis.ad}
\alias{gw.weight.box}
\alias{gw.weight.box.ad}
\alias{gw.weight.gau}
\alias{gw.weight.gau.ad}
\alias{gw.weight.tri}
\alias{gw.weight.tri.ad}
\alias{bisq_wt_vec}
\alias{bisq_wt_mat}
\alias{gauss_wt_vec}
\alias{gauss_wt_mat}
\alias{tri_wt_vec}
\alias{tri_wt_mat}
\alias{exp_wt_vec}
\alias{exp_wt_mat}
\title{Weight matrix calculation}
\description{
  Calculate a weight vector(matrix) from a distance vector(matrix).
}
\usage{
gw.weight(vdist,bw,kernel,adaptive=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vdist}{a distance matrix or vector}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
}

\value{
 Returns a numeric weight matrix or vector;
 matrix with its rows corresponding to the observations and its columns corresponds to the GW model calibration points.
}
\note{
The gaussian and exponential kernel functions are continuous and valued in the interval (0,1]; while bisquare, tricube and 
boxcar kernel functions are discontinuous and valued in the interval [0,1]. Notably, the upper limit of the bandwidth is exactly
the number of observations when the adaptive kernel is used. In this function, the adaptive bandwidth will be specified as the number
of observations even though a larger number is assigned. The function will be the same as a global application function (i.e. all weights are 1) when 
the adaptive bandwidth is equal to or larger than the number of observations when using the boxcar kernel function.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GW tools}
\concept{weight}
