\name{Georgia}
\alias{Georgia}
\alias{Gedu.df}
\docType{data}
\title{Georgia census data set (csv file)}
\description{
  Census data from the county of Georgia, USA 
}
\usage{data(Georgia)}
\format{
  A data frame with 159 observations on the following 13 variables.
  \describe{
    \item{AreaKey}{An identification number for each county}
    \item{Latitude}{The latitude of the county centroid}
    \item{Longitud}{The longitude of the county centroid}
    \item{TotPop90}{Population of the county in 1990}
    \item{PctRural}{Percentage of the county population defined as rural}
    \item{PctBach}{Percentage of the county population with a bachelors degree}
    \item{PctEld}{Percentage of the county population aged 65 or over}
    \item{PctFB}{Percentage of the county population born outside the US}
    \item{PctPov}{Percentage of the county population living below the poverty line}
    \item{PctBlack}{Percentage of the county population who are black}
    \item{ID}{a numeric vector of IDs}
    \item{X}{a numeric vector of x coordinates}
    \item{Y}{a numeric vector of y coordinates}
  }
}
\details{
This data set can also be found in GWR 3 and in spgwr.
}
\references{
  Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.
}
\examples{
data(Georgia)
ls()
coords <- cbind(Gedu.df$X, Gedu.df$Y)
educ.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
spplot(educ.spdf, names(educ.spdf)[4:10])
}
\keyword{data}
\concept{Georgia census}
