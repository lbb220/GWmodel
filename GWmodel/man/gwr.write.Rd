\name{gwr.write}
\alias{gwr.write}
\alias{writeGWR}
\alias{writeGWR.shp}
\alias{gwr.write.shp}
\title{Write the GWR results into files}
\description{
This function writes the calibration result of function \link{gwr.basic} to a text file and shape files}
\usage{
gwr.write(x,fn="GWRresults")
gwr.write.shp(x,fn="GWRresults")
}
\arguments{
  \item{x}{an object of class \dQuote{gwrm}, returned by the function \link{gwr.basic}}
  \item{fn}{file name for the written results, by default the output files can be found in the working directory, \dQuote{GWRresults.txt}, \dQuote{GWRresults(.shp, .shx, .dbf)}}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\note{
  The projection file is missing for the writen shapefiles.
  
  The functions \dQuote{writeGWR} and \dQuote{writeGWR.shp} (in the early versions of GWmodel) have been renamed respectively as
 \dQuote{gwr.write} and \dQuote{gwr.write.shp}, while the old names are still kept valid.
}
\keyword{GWR}
\concept{results writing}