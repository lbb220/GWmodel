\name{gwr.model.view}
\alias{gwr.model.view}
\alias{model.view.gwr}
\title{Visualise the GWR models from \code{\link{gwr.model.selection}}}
\description{
This function visualises the GWR models from \code{\link{gwr.model.selection}}.
}
\usage{
gwr.model.view(DeVar, InDeVars, model.list)
}
\arguments{
  \item{DeVar}{dependent variable}
  \item{InDeVars}{a vector of independent variables for model selection}
  \item{model.list}{a list of all GWR model tried in \code{\link{gwr.model.selection}}}
}

\note{
  The function \dQuote{model.view.gwr} (in the early versions of GWmodel) has been renamed as
 \dQuote{gwr.model.view}, while the old name is still kept valid.
}

\author{Binbin Lu \email{binbinlu@whu.edu.cn}}

\seealso{\code{\link{gwr.model.selection}}, \code{\link{gwr.model.sort}}}
\examples{
\dontrun{
data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
DeVar<-"PURCHASE"
InDeVars<-c("FLOORSZ","GARAGE1","BLDPWW1","BLDPOSTW")
model.sel<-gwr.model.selection(DeVar,InDeVars, data=londonhp,
kernel = "gaussian", dMat=DM,bw=5000)
model.list<-model.sel[[1]]
gwr.model.view(DeVar, InDeVars, model.list=model.list)
}
}
\keyword{GWR}
\concept{model selection}

