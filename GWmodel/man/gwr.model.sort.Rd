\name{gwr.model.sort}
\alias{gwr.model.sort}
\alias{model.sort.gwr}
\title{Sort the results of the GWR model selection function \code{\link{gwr.model.selection}}.}
\description{
Sort the results from the GWR model selection function \code{\link{gwr.model.selection}}
}
\usage{
gwr.model.sort(Sorting.list , numVars, ruler.vector)
}
\arguments{
  \item{Sorting.list}{a list returned by function \code{\link{gwr.model.selection}}}
  \item{numVars}{the number of independent variables involved in model selection}
  \item{ruler.vector}{a numeric vector as the sorting basis}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\note{
  The function sorts the results of model selection within individual levels.
  
  The function \dQuote{model.sort.gwr} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwr.model.sort}, while the old name is still kept valid.
}
\seealso{\code{\link{gwr.model.selection}}, \code{\link{gwr.model.view}}}
\keyword{GWR}
\concept{model selection}

