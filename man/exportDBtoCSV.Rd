\name{exportDBtoCSV}
\alias{exportDBtoCSV}
\title{
Export a LOBSTAHS database (LOBdbase) object to a text file
}
\description{
Collects all elements of a \code{\link{LOBdbase-class}} object into a single
table, then writes the table to a .csv file.
}
\usage{
exportDBtoCSV(LOBdbase)
}
\arguments{
  \item{LOBdbase}{
A LOBSTAHS \code{"LOBdbase"} object.
}
}
\details{\code{exportDBtoCSV} extracts data from all slots in a given 
\code{"LOBdbase"} object and flows it into a data frame. The data frame is then
written to a .csv file with a series of plain-language column headers. The file
name includes the database polarity mode and a timestamp (in user's local time)
indicating the time and date of export.
}
\value{
  A .csv file containing data from a \code{\link{LOBdbase-class}} object.
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy. 2016.
LOBSTAHS: An adduct-based lipidomics strategy for discovery and identification
of oxidative stress biomarkers. \emph{Analytical Chemistry} 88:7154-7162 
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\link{LOBdbase}},
  \code{\linkS4class{LOBdbase}},
  \code{\link{generateLOBdbase}}

}
\examples{
## export default positive mode database to a text file

data("default.LOBdbase")

exportDBtoCSV(default.LOBdbase$positive)

## file is written to user's working directory, with filename containing
## polarity and timestamp of generation
}
\keyword{methods}
