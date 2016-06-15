\name{loadLOBdbase}
\alias{loadLOBdbase}
\title{
Import and reconstruct LOBdbase from a text file
}
\description{
Reconstruct a \code{"LOBdbase"} object from properly formatted data in a .csv
file. More robust than the generic constructor function
\code{\link{LOBdbase}}.
}
\usage{
loadLOBdbase(file, polarity, num.compounds = NULL)
}
\arguments{
  \item{file}{
Path to a .csv file containing the database to be imported. The file format
should be consistent with the .csv output that is obtained using
\code{gen.csv = TRUE} in \code{\link{generateLOBdbase}}. Some leeway is allowed
for variation in capitalization and punctuation of column headers.
}
  \item{polarity}{
Specify polarity (ionization mode) of the database being imported ("positive" or
"negative"). \code{loadLOBdbase} will attempt to detect the polarity mode based
on data in the "adduct" column if nothing is given.
}
  \item{num.compounds}{
If known, the number of parent compounds represented in the database being
imported. Can be unspecified.
}
}
\details{\code{loadLOBdbase} is a glorified implementation of
\code{\link[utils]{read.table}} that attempts to determine whether the format of
the data in \code{file} is consistent with that of a saved \code{"LOBdbase"}
object. Some leeway is allowed for variation in punctuation and spelling of
column headers in the source file. Special attention should be paid to how
adduct ions are specified (e.g., \code{"[M+H]+"}, \code{"[M+NH4]+"}, or
\code{"[M+Cl]-"}). Most users will find this function more useful than the
generic constructor function \code{\link{LOBdbase}}.
}
\value{
  A \code{"\link{LOBdbase-class}"} object.
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy. 2016.
LOBSTAHS: An adduct-based lipidomics strategy for discovery and identification
of oxidative stress biomarkers. Analytical Chemistry.
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\link{LOBdbase}},
  \code{\linkS4class{LOBdbase-class}},
  \code{\link{doLOBscreen}},
  \code{\link{generateLOBdbase}},
  \code{\link{default.LOBdbase}}
}
\examples{
## save the default negative mode database as a .csv file

data(default.LOBdbase)

neg.DB = default.LOBdbase$negative

fname = paste0("LOBSTAHS_lipid-oxy_DB_",
               strtrim(as.character(polarity(neg.DB)),3),".csv")

exportmat = data.frame(frag_ID(neg.DB),
                         mz(neg.DB),
                         exact_parent_neutral_mass(neg.DB),
                         as.character(lipid_class(neg.DB)),
                         as.character(species(neg.DB)),
                         as.character(adduct(neg.DB)),
                         as.character(adduct_rank(neg.DB)),
                         FA_total_no_C(neg.DB),
                         FA_total_no_DB(neg.DB),
                         degree_oxidation(neg.DB),
                         parent_elem_formula(neg.DB),
                         parent_compound_name(neg.DB),
                         stringsAsFactors = FALSE)

colnames(exportmat) = c("frag_ID","mz","exact_parent_neutral_mass",
                        "lipid_class","species","adduct","adduct_rank",
                        "FA_total_no_C","FA_total_no_DB","degree_oxidation",
                        "parent_elem_formula","parent_compound_name")

write.csv(exportmat, fname)

## reimport it

neg.DB.reimported = loadLOBdbase("LOBSTAHS_lipid-oxy_DB_neg.csv",
                                 polarity = "negative",
                                 num.compounds = NULL)
}
\keyword{methods}
