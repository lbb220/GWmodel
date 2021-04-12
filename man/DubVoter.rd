\name{DubVoter}
\alias{DubVoter}
\alias{Dub.voter}
\docType{data}
\title{Voter turnout data in Greater Dublin(SpatialPolygonsDataFrame)}
\description{
  Voter turnout and social characters data in Greater Dublin for the 2002 General election and the 2002 census.
  Note that this data set was originally thought to relate to 2004, so for continuity we have retained the associated variable names.
}
\usage{data(DubVoter)}
\format{
  A SpatialPolygonsDataFrame with 322 electoral divisions on the following 11 variables.
  \describe{
    \item{DED_ID}{a vector of ID}
    \item{X}{a numeric vector of x coordinates}
    \item{Y}{a numeric vector of y coordinates}
    \item{DiffAdd}{percentage of the population in each ED who are one-year migrants (i.e. moved to a different address 1 year ago)}
    \item{LARent}{percentage of the population in each ED who are local authority renters}
    \item{SC1}{percentage of the population in each ED who are social class one (high social class)}
    \item{Unempl}{percentage of the population in each ED who are unemployed}
    \item{LowEduc}{percentage of the population in each ED who are with little formal education}
    \item{Age18_24}{percentage of the population in each ED who are age group 18-24}
    \item{Age25_44}{percentage of the population in each ED who are age group 25-44}
    \item{Age45_64}{percentage of the population in each ED who are age group 45-64}
    \item{GenEl2004}{percentage of population in each ED who voted in 2004 election}    
  }
}
\details{
Variables are from DubVoter.shp.
}
\references{
  Kavanagh A (2006) Turnout or turned off? Electoral participation in Dublin in the early 21st Century. Journal of Irish Urban Studies 3(2):1-24
  
  Harris P, Brunsdon C, Charlton M (2011) Geographically weighted principal components analysis.  International Journal of Geographical Information Science 25 (10):1717-1736
}
\examples{

data(DubVoter)
ls()
\dontrun{
spplot(Dub.voter,names(Dub.voter)[4:12])
}
}
\keyword{data}
\concept{Dublin Voter turnout}
