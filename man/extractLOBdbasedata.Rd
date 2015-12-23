\name{extractLOBdbasedata}
\alias{extractLOBdbasedata}
\title{
Extract Compound Data from a LOBdbase
}
\description{
Extracts data from a LOBSTAHS \code{"\link{LOBdbase}"} object for a given molecular species.
}
\usage{
extractLOBdbasedata(frag_ID, database)
}
\arguments{
  \item{frag_ID}{
The LOBdbase fragment ID(s) corresponding to the species for which data should be extracted.
}
  \item{database}{
The \code{LOBdbase} from which data are to be extracted.
}
}
\value{
A data frame containing the extracted data, with field names corresponding to the slots in \code{database}.
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\linkS4class{LOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{generateLOBdbase}},
  \code{\link{default.LOBdbase}}  
}
\examples{
## extract data from the default database for the [M+H]+ adduct of
## dinitrophenyl-phosphatidylethanolamine (DNPPE), a synthetic lipid used as an
## internal standard

library(LOBSTAHS)

data(default.LOBdbase)

DNPPE.data = extractLOBdbasedata(frag_ID = 19801, database = LOBdbase$positive)

}
\keyword{methods}
