\name{GWmodel-package}
\alias{GWmodel-package}
\alias{GWmodel}
\docType{package}
\title{
Geographically-Weighted Models
}
\description{
In GWmodel, we introduce techniques from a particular branch of spatial statistics, 
termed geographically-weighted (GW) models. GW models suit situations when data are 
not described well by some global model, but where there are spatial regions where 
a suitably localised calibration provides a better description. GWmodel includes functions 
to calibrate: GW summary statistics, GW principal components analysis, GW discriminant analysis 
and various forms of GW regression; some of which are provided in basic and robust (outlier resistant) forms.
In particular, the high-performence computing technologies, including multi-thread and CUDA techniques are started
to be adopted for efficient calibrations.
}
\details{
\tabular{ll}{
Package: \tab GWmodel\cr
Type: \tab Package\cr
Version: \tab 2.2-4\cr
Date: \tab 2021-02-13\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Binbin Lu, Paul Harris, Martin Charlton, Chris Brunsdon, Tomoki Nakaya, Daisuke Murakami,Isabella Gollini[ctb], Yigong Hu[ctb], Fiona H Evans[ctb]

Maintainer: Binbin Lu <binbinlu@whu.edu.cn>
}
\note{
Acknowledgements: We gratefully acknowledge support from National Natural Science Foundation of China (42071368); Science Foundation Ireland 
under the National Development Plan through the award of a Strategic Research Centre 
grant 07-SRC-I1168.

Beta versions can always be found at \url{https://github.com/lbb220/GWmodel},
which includes all the newly developed functions for GW models.

For latest tutorials on using GWmodel please go to: \url{https://rpubs.com/gwmodel}
}

\references{
Gollini I, Lu B, Charlton M, Brunsdon C, Harris P (2015) GWmodel: an R Package 
for exploring Spatial Heterogeneity using Geographically Weighted Models. Journal 
of Statistical Software, 63(17):1-50, \url{https://www.jstatsoft.org/article/view/v063i17}

Lu B, Harris P, Charlton M, Brunsdon C (2014) The GWmodel R Package: further 
topics for exploring Spatial Heterogeneity using Geographically Weighted Models.
Geo-spatial Information Science 17(2): 85-101, 
\url{https://www.tandfonline.com/doi/abs/10.1080/10095020.2014.917453}
}
\keyword{package}
