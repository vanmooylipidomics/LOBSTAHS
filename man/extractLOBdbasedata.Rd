\name{extractLOBdbasedata}
\alias{extractLOBdbasedata}
\title{
Extract compound data from a LOBdbase
}
\description{
Extracts data from a LOBSTAHS \code{\link{LOBdbase-class}} object for a given 
molecular species.
}
\usage{
extractLOBdbasedata(frag_ID, database)
}
\arguments{
  \item{frag_ID}{
The LOBdbase fragment ID(s) corresponding to the species for which data should 
be extracted.
}
  \item{database}{
The \code{LOBdbase} from which data are to be extracted.
}
}
\value{
A data frame containing the extracted data, with field names corresponding to 
the slots in \code{database}.
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
  \code{\linkS4class{LOBdbase}},
  \code{\link{LOBdbase}}
  \code{\link{loadLOBdbase}},
  \code{\link{generateLOBdbase}},
  \code{\link{default.LOBdbase}}  
}
\examples{
## extract data from the default database for the [M+H]+ adduct of
## dinitrophenyl-phosphatidylethanolamine (DNPPE), a synthetic lipid used as an
## internal standard

data(default.LOBdbase)

DNPPE.data = extractLOBdbasedata(frag_ID = 19801,
                                 database = default.LOBdbase$positive)

}
\keyword{methods}
