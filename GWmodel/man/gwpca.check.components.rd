\name{gwpca.check.components}
\alias{gwpca.check.components}
\alias{check.components}
\title{Interaction tool with the GWPCA glyph map}
\description{
The function interacts with the multivariate glyph plot of GWPCA loadings.
}
\usage{
gwpca.check.components(ld,loc)
}

\arguments{
  \item{ld}{GWPCA loadings returned by \link{gwpca}}
  \item{loc}{a 2-column numeric array of GWPCA evaluation locations}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\note{
The function \dQuote{check.components} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwpca.check.components}, while the old name is still kept valid.
}
\seealso{\link{gwpca.glyph.plot}}
\keyword{GWPCA}
\concept{glyph plot interaction}
