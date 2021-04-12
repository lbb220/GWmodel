\name{GE2015}
\alias{ge2015}
\docType{data}
\title{General Election outcome data at constituency level in England}
\description{
  2015 General Election outcome data and associated socio-economic indicators at constituency level in England
}
\usage{data(GE2015)}
\format{
  A SpatialPolygonsDataFrame with the following variables of interest.
  \describe{
    \item{WINNER}{Elected party in each constituency in 2015}
    \item{Winner10}{Elected party in each constituency in 2010}
    \item{Winner15}{elected party in each constituency in 2015}
    \item{Age65over}{Percentage of the population aged 65 or over}
    \item{OwnOcc}{Percentage of households either fully owned or mortgaged by the occupier}
    \item{NoQual}{Percentage of population with no educational qualifications} 
    \item{NonWhite}{Percentage of the population who are not white} 
    \item{LoneParHH}{Percentage of households whose head is a lone parent}
    \item{Unemp}{Percentage of households whose designated head of house is economically active but not employed.} 	
  }
}
\examples{
data(GE2015)
ls()
}
\keyword{data}
\concept{General Election outcome data in England}
