\name{gwpca.glyph.plot}
\alias{gwpca.glyph.plot}
\alias{glyph.plot}
\title{Multivariate glyph plots of GWPCA loadings}
\description{
This function provides a multivariate glyph plot of GWPCA loadings at each output location.
}
\usage{
gwpca.glyph.plot(ld,loc, r1=50, add=FALSE,alpha=1,sep.contrasts=FALSE)
}

\arguments{
  \item{ld}{GWPCA loadings returned by \link{gwpca}}
  \item{loc}{a two-column numeric array for providing evaluation locations of GWPCA calibration}
  \item{r1}{argument for the size of the glyphs, default is 50; glyphs get larger as r1 is reduced}
  \item{add}{if TRUE, add the plot to the existing window.}
  \item{alpha}{the level of transparency of glyph from function rgb() and ranges from 0 to max (fully transparent to opaque)}
  \item{sep.contrasts}{allows different types of glyphs and relates to whether absolute loadings are used (TRUE) or not}
}
\note{
The function \dQuote{glyph.plot} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwpca.glyph.plot}, while the old name is still kept valid.
}
\references{
Harris P, Brunsdon C, Charlton M (2011)
Geographically weighted principal components analysis.
International Journal of Geographical Information Science 25:1717-1736
}
\keyword{GWPCA}
\concept{glyph plot}

