\name{gwr.mink.matrixview}
\alias{gwr.mink.matrixview}
\alias{mink.matrixview}
\title{Visualisation of the results from \code{\link{gwr.mink.approach}}}
\description{
This function visualises the AICc/CV results from the \code{\link{gwr.mink.approach}}.
}
\usage{
gwr.mink.matrixview(diag.df, znm=colnames(diag.df)[4], criterion="AIC")
}

\arguments{
  \item{diag.df}{the first part of a list object returned by \code{\link{gwr.mink.approach}}}
  \item{znm}{the name of the forth column in diag.df}
  \item{criterion}{the criterion used for distance metric selection in \code{\link{gwr.mink.approach}}}
}
\note{
The function \dQuote{mink.matrixview} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwr.mink.matrixview}, while the old name is still kept valid.
}
\references{
Lu, B, Charlton, M, Brunsdon, C & Harris, P(2016). The Minkowski approach for choosing the distance metric in Geographically Weighted Regression. International Journal of Geographical Information Science, 30(2): 351-368.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{GWR}
\concept{Minkowski approach view}

