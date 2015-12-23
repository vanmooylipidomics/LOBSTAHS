\name{extractLOBdbasedata}
\alias{extractLOBdbasedata}
\title{
Extract Compound Data from a LOBdbase
}
\description{
Extracts data from a LOBSTAHS \code{\link[LOBSTAHS]{"LOBdbase"}} object for a given compound or compounds.
}
\usage{
extractLOBdbasedata(frag_ID, database)
}
\arguments{
  \item{frag_ID}{
The LOBdbase fragment ID(s) corresponding to the compounds for which data should be extracted.
}
  \item{database}{
The \code{LOBdbase} from which data are to be extracted.
}
}
\value{
A data frame containing the extracted data, with field names corresponding to the slots in \code{database}.
}
\references{
\item{Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\linkS4class[LOBSTAHS]{LOBdbase}},
  \code{\link[LOBSTAHS]{loadLOBdbase}},
  \code{\link[LOBSTAHS]{generateLOBdbase}},
  \code{\link[LOBSTAHS]{default.LOBdbase}}  
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
