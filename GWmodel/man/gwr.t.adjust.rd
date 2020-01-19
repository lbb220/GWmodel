\name{gwr.t.adjust}
\alias{gwr.t.adjust}
\title{Adjust p-values for multiple hypothesis tests in basic GWR}
\description{
Given a set of p-values from the pseudo t-tests of basic GWR outputs, this function returns adjusted p-values using: 
(a) Bonferroni, (b) Benjamini-Hochberg, (c) Benjamini-Yekutieli and
(d) Fotheringham-Byrne procedures.
}
\usage{
gwr.t.adjust(gwm.Obj)
}
\arguments{
  \item{gwm.Obj}{an object of class \dQuote{gwrm}, returned by the function \link{gwr.basic}}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\references{
Byrne, G., Charlton, M. and Fotheringham, S., 2009. Multiple dependent hypothesis tests 
in geographically weighted regression. In: Lees, B. and Laffan, S. eds. 10th 
International conference on geocomputation. Sydney.
}
\keyword{GWR}
\keyword{p-values adjustment}
