\name{LOBdbase}
\Rdversion{1.2}
\alias{LOBdbase}

\title{LOBdbase constructor for manual creation or reconstruction of a LOBdbase
object}
\description{
Constructor function for manual assembly or reconstruction of a
\code{LOBdbase} object to be used in screening data with the
\code{\link{LOBSTAHS}} function \code{\link{doLOBscreen}}.
}
\usage{
LOBdbase(frag_ID = NULL, mz = NULL, exact_parent_neutral_mass = NULL,
        lipid_class = NULL, species = NULL, adduct = NULL, adduct_rank = NULL,
        FA_total_no_C = NULL, FA_total_no_DB = NULL, degree_oxidation = NULL,
        parent_elem_formula = NULL, parent_compound_name = NULL,
        polarity = NULL, num_entries = NULL, num_compounds = NULL)
}
\arguments{
  \item{frag_ID}{
An object of class \code{"integer"}; vector of unique identifiers for the
molecular species in the database.
}
  \item{mz}{
An object of class \code{"numeric"}; calculated \emph{m/z} of each species for
which an entry exists.
}
  \item{exact_parent_neutral_mass}{
An object of class \code{"numeric"}; calculated (monoisotopic) exact masses of
the parent compound of each species.
}
   \item{exact_parent_neutral_mass}{
An object of class \code{"factor"}; parent lipid classes of each species
}
   \item{species}{
An object of class \code{"character"}; the lipid subclasses of each species
}
   \item{adduct}{
An object of class \code{"factor"}; the adduct ions represented by each entry
}
   \item{adduct_rank}{
An object of class \code{"integer"}; the relative abundance rankings of each
adduct relative to the other adducts of the same parent compound
}
   \item{FA_total_no_C}{
An object of class \code{"integer"}; the total number  of acyl (fatty acid)
carbon atoms in the parent compound of each entry; values should be \code{NA}
where \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA
}
   \item{FA_total_no_DB}{
An object of class \code{"integer"}; the total number of acyl (fatty acid)
carbon-carbon double bonds in the parent compound of each entry; values should
be \code{NA} where \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA
}
   \item{degree_oxidation}{
An object of class \code{"integer"}; the number of additional oxygen atoms
present on each species
}
   \item{parent_elem_formula}{
An object of class \code{"character"}; the elemental formulae of the parent
compound of each species
}
   \item{parent_compound_name}{
An object of class \code{"character"}; names of the parent compound of each
species; see the reference for this entry for snaming conventions that should
be applied for pigments and compounds other than pigments
}
   \item{polarity}{
An object of class \code{"factor"} with length = 1; the ionization mode of data
in the database to be constructed
}
   \item{num_entries}{
An object of class \code{"integer"} with length = 1; the number of total entries
(adducts) in the database to be constructed. If specified, the value of
\code{num_entries} should be equal to the length of any arguments that contain
database data
}
   \item{num_compounds}{
An object of class \code{"integer"} with length = 1; the number of parent
compounds represented in the database to be constructed. If specified, the
value of \code{num_compounds} should be < the value specifed for
\code{num_entries}.
}
}
\details{Typically, a \code{LOBdbase} will be created using the in silico
simulation function \code{\link{generateLOBdbase}}. Formatted database entries
(such as from an external .csv file) can be loaded using
\code{\link{loadLOBdbase}}. The rudimentary  \code{LOBdbase} constructor
function is therefore provided only for manual object creation; it will not be
needed by most users. All arguments except for those containing metadata
(i.e., \code{polarity}, \code{num_entries}, and \code{num_compounds}) should be
of the same length.
}
\value{
  A \code{"LOBdbase-class"} object.
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
  \code{\link{generateLOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{doLOBscreen}},
  \code{\link{LOBdbase}}
}
\examples{
## create an empty LOBdbase

myLOBdbase = LOBdbase(frag_ID = NULL, mz = NULL,
                      exact_parent_neutral_mass = NULL, lipid_class = NULL,
                      species = NULL, adduct = NULL, adduct_rank = NULL,
                      FA_total_no_C = NULL, FA_total_no_DB = NULL,
                      degree_oxidation = NULL, parent_elem_formula = NULL,
                      parent_compound_name = NULL, polarity = NULL,
                      num_entries = NULL, num_compounds = NULL)
}
\keyword{methods}
