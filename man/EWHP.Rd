\name{EWHP}
\alias{EWHP}
\alias{ewhp}
\docType{data}
\title{House price data set (DataFrame) in England and Wales}
\description{
A house price data set for England and Wales from 2001 with 9 hedonic (explanatory) variables.
}
\usage{data(EWHP)}
\format{
 A data frame with 519 observations on the following 12 variables.
  \describe{
    \item{Easting}{a numeric vector, X coordinate}
    \item{Northing}{a numeric vector, Y coordinate}
    \item{PurPrice}{a numeric vector, the purchase price of the property}
    \item{BldIntWr}{a numeric vector, 1 if the property was built during the world war, 0 otherwise}
    \item{BldPostW}{a numeric vector, 1 if the property was built after the world war, 0 otherwise}
    \item{Bld60s}{a numeric vector, 1 if the property was built between 1960 and 1969, 0 otherwise}
    \item{Bld70s}{a numeric vector, 1 if the property was built between 1970 and 1979, 0 otherwise}
    \item{Bld80s}{a numeric vector, 1 if the property was built between 1980 and 1989, 0 otherwise}
    \item{TypDetch}{a numeric vector, 1 if the property is detached (i.e. it is a stand-alone house), 0 otherwise}
    \item{TypSemiD}{a numeric vector, 1 if the property is semi detached, 0 otherwise}
    \item{TypFlat}{a numeric vector, if the property is a flat (or 'apartment' in the USA), 0 otherwise}
    \item{FlrArea}{a numeric vector, floor area of the property in square metres}
  }
}
\references{
  Fotheringham, A.S., Brunsdon, C., and Charlton, M.E. (2002), Geographically Weighted Regression: 
  The Analysis of Spatially Varying Relationships, Chichester: Wiley.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
###
data(EWHP)
head(ewhp)
houses.spdf <- SpatialPointsDataFrame(ewhp[, 1:2], ewhp)
 ####Get the border of England and Wales
data(EWOutline)
plot(ewoutline)
plot(houses.spdf, add = TRUE, pch = 16)
}
\keyword{data}
\concept{house price}
