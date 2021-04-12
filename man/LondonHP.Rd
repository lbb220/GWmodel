\name{LondonHP}
\alias{LondonHP}
\alias{londonhp}
\docType{data}
\title{London house price data set (SpatialPointsDataFrame)}
\description{
A house price data set with 18 hedonic variables for London in 2001.
}
\usage{data(LondonHP)}
\format{
 A SpatialPointsDataFrame object (proj4string set to "+init=epsg:27700 +datum=OSGB36").

  The "data" slot is a data frame with 372 observations on the following 21 variables.
  \describe{
    \item{X}{a numeric vector, X coordinate}
    \item{Y}{a numeric vector, Y coordinate}
    \item{PURCHASE}{a numeric vector, the purchase price of the property}
    \item{FLOORSZ}{a numeric vector, floor area of the property in square metres}
    \item{TYPEDETCH}{a numeric vector, 1 if the property is detached (i.e. it is a stand-alone house), 0 otherwise}
    \item{TPSEMIDTCH}{a numeric vector, 1 if the property is semi detached, 0 otherwise}
    \item{TYPETRRD}{a numeric vector, 1 if the property is in a terrace of similar houses (commonly referred to as a 'row house' in the USA), 0 otherwise}
    \item{TYPEBNGLW}{a numeric vector, if the property is a bungalow (i.e. it has only one floor), 0 otherwise}
    \item{TYPEFLAT}{a numeric vector, if the property is a flat (or 'apartment' in the USA), 0 otherwise}
    \item{BLDPWW1}{a numeric vector, 1 if the property was built prior to 1914, 0 otherwise}
    \item{BLDPOSTW}{a numeric vector, 1 if the property was built between 1940 and 1959, 0 otherwise}
    \item{BLD60S}{a numeric vector, 1 if the property was built between 1960 and 1969, 0 otherwise}
    \item{BLD70S}{a numeric vector, 1 if the property was built between 1970 and 1979, 0 otherwise}
    \item{BLD80S}{a numeric vector, 1 if the property was built between 1980 and 1989, 0 otherwise}
    \item{BLD90S}{a numeric vector, 1 if the property was built between 1990 and 2000, 0 otherwise}
    \item{BATH2}{a numeric vector, 1 if the property has more than 2 bathrooms, 0 otherwise}
    \item{GARAGE}{a numeric vector,1 if the house has a garage, 0 otherwise}
    \item{CENTHEAT}{a numeric vector, 1 if the house has central heating, 0 otherwise}
    \item{BEDS2}{a numeric vector, 1 if the property has more than 2 bedrooms, 0 otherwise}
    \item{UNEMPLOY}{a numeric vector, the rate of unemployment in the census ward in which the house is located}
    \item{PROF}{a numeric vector, the proportion of the workforce in professional or managerial occupations in the census ward in which the house is located}
  }
}
\references{
Fotheringham, A.S., Brunsdon, C., and Charlton, M.E. (2002), Geographically Weighted Regression: 
The Analysis of Spatially Varying Relationships, Chichester: Wiley.
  
Lu, B, Charlton, M, Harris, P, Fotheringham, AS (2014) Geographically weighted regression 
with a non-Euclidean distance metric: a case study using hedonic house price data. 
International Journal of Geographical Information Science 28(4): 660-681
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\examples{
data(LondonHP)
data(LondonBorough)
ls()
plot(londonborough)
plot(londonhp, add=TRUE)
}
\keyword{data}
\concept{house price}
