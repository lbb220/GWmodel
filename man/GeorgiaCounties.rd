\name{GeorgiaCounties}
\alias{Gedu.counties}
\docType{data}
\title{Georgia counties data (SpatialPolygonsDataFrame)}
\description{
  The Georgia census data with boundaries for mapping
}
\usage{data(GeorgiaCounties)}
\details{
This data set can also be found in GWR 3 and in spgwr.
}
\examples{
data(GeorgiaCounties)
plot(Gedu.counties)
data(Georgia)
coords <- cbind(Gedu.df$X, Gedu.df$Y)
educ.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
plot(educ.spdf, add=TRUE)

}
\keyword{data}
\concept{Georgia counties}
