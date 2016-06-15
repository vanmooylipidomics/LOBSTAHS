\name{generateLOBdbase}
\alias{generateLOBdbase}
\title{
Conduct \emph{in silico} simulation and generate lipid-oxylipin database
}
\description{
Applies an \emph{in silico} simulation to generate data by ionization mode
(polarity) for a wide range of lipids, oxidized lipids, and oxylipins. 
User-supplied structural criteria and empirically-determined adduct ion 
abundance rankings for the major lipid classes are used to create entries for a 
range of lipid moieties. The database(s) can then be used in 
\code{\link{doLOBscreen}} to assign compound identities to grouped peakdata.
}
\usage{
generateLOBdbase(polarity = c("positive","negative"), gen.csv = FALSE,
                 component.defs = NULL, AIH.defs = NULL, acyl.ranges = NULL,
                 oxy.ranges = NULL)
}
\arguments{
  \item{polarity}{
Ionization mode for which database is to be generated.
}
  \item{gen.csv}{
Should results also be written to a .csv file?
}
  \item{component.defs}{
File path to a .csv file containing elemental composition definitions for the 
various chemical components needed by \code{generateLOBdbase}. If nothing is 
specified, \code{generateLOBdbase} will use the default composition table 
(\code{\link{default.componentCompTable}}). The default table includes 
definitions for the masses of a wide range of adducts, photosynthetic pigments, 
and structural backbones of some major lipid classes.

A Microsoft Excel spreadsheet template included with the package at 
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_basic_component_matrix.xlsx can be 
used to generate a .custom csv file with elemental composition definitions in a 
format appropriate for \code{generateLOBdbase}. Alternatively, the spreadsheet 
may be downloaded at \url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_basic_component_matrix.xlsx}.
}
  \item{AIH.defs}{
File path to a .csv file containing empirical adduct ion hierarchy (AIH) data 
for various pigments and lipid classes. If nothing is specified, 
\code{generateLOBdbase} will use the default AIH data in 
(\code{\link{default.adductHierarchies}}).

A Microsoft Excel spreadsheet template included with the package at 
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_adduct_ion_hierarchies.xlsx can be 
used to generate a .csv file with additional (or alternative) adduct hierarchy 
data in a format appropriate for \code{generateLOBdbase}. Alternatively, the 
spreadsheet may be downloaded at \url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_adduct_ion_hierarchies.xlsx}.
}
  \item{acyl.ranges}{
File path to a .csv file containing ranges of values for the total number of 
acyl (i.e., fatty acid) carbon atoms to be considered during the 
\emph{in silico} similation on intact polar diacylglycerols (IP-DAG), 
triacylglycerols (TAG), polyunsaturated aldehydes (PUAs), and free fatty acids 
(FFA). If nothing is specified, \code{generateLOBdbase} will use the default 
acyl carbon atom range data in (\code{\link{default.acylRanges}}).

A Microsoft Excel spreadsheet template included with the package at 
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_acyl_prop_ranges.xlsx can be used 
to generate a .csv file with custom acyl carbon range data in a format 
appropriate for \code{generateLOBdbase}. Alternatively, the spreadsheet may be 
downloaded at \url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_acyl_prop_ranges.xlsx}.
}
  \item{oxy.ranges}{
File path to a .csv file containing ranges of values for the number of 
additional oxygen atoms to be considered during the \emph{in silico} similation 
on intact polar diacylglycerols (IP-DAG), triacylglycerols (TAG), 
polyunsaturated aldehydes (PUAs), and free fatty acids (FFA). If nothing is 
specified, \code{generateLOBdbase} will use the default oxidation state ranges 
in (\code{\link{default.oxyRanges}}).

A Microsoft Excel spreadsheet template included with the package at 
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_addl_oxy_ranges.xlsx can be used to
generate a .csv file with custom oxidation state ranges in a format appropriate 
for \code{generateLOBdbase}. Alternatively, the spreadsheet may be downloaded at
\url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_addl_oxy_ranges.xlsx}. By default, \code{generateLOBdbase} considers 0-4 additional oxygen 
atoms on each chemically possible IP-DAG, TAG, PUA, and FFA.
}
}
\details{Using the default structural property inputs described here, 
\code{generateLOBdbase} can produce databases with entries for a wide range of 
intact polar diacylglycerols (IP-DAG), triacylglycerols (TAG), polyunsaturated 
aldehydes (PUAs), free fatty acids (FFA), and common photosynthetic pigments. 
The default databases (as of June 2016) contain data on 13,578 and 11,073 unique
compounds that can be identifed in positive and negative ion mode spectra, 
respectively.

\strong{Note} that the default databases have been pre-generated (see 
\code{\link{default.LOBdbase}}) and it is therefore unnecessary to call 
\code{generateLOBdbase} with the default parameters.
}
\value{
  A \code{"\link{LOBdbase}"} object with the structure:
  \describe{
    \item{\code{frag_ID}:}{Object of class \code{"integer"}, a unique identifier
    for this molecular species }
    \item{\code{mz}:}{Object of class \code{"numeric"}, the calculated 
    \emph{m/z} of this species }
    \item{\code{exact_parent_neutral_mass}:}{Object of class \code{"numeric"}, 
    the calculated (monoisotopic) exact mass of the parent compound of this 
    species }
    \item{\code{lipid_class}:}{Object of class \code{"factor"}, the parent lipid
    class of this species }
    \item{\code{species}:}{Object of class \code{"character"}, the lipid 
    subclass }
    \item{\code{adduct}:}{Object of class \code{"factor"}, the adduct ion 
    represented by this entry }
    \item{\code{adduct_rank}:}{Object of class \code{"integer"}, the relative 
    abundance ranking of this adduct relative to the other adducts of the same 
    parent compound }
    \item{\code{FA_total_no_C}:}{Object of class \code{"integer"}, total number 
    of acyl (fatty acid) carbon atoms in the parent compound; \code{NA} if 
    \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{FA_total_no_DB}:}{Object of class \code{"integer"}, total number
    of acyl (fatty acid) carbon-carbon double bonds in the parent compound; 
    \code{NA} if \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{degree_oxidation}:}{Object of class \code{"integer"}, number of 
    additional oxygen atoms present }
    \item{\code{parent_elem_formula}:}{Object of class \code{"character"}, 
    elemental formula of the parent compound }
    \item{\code{parent_compound_name}:}{Object of class \code{"character"}, name
    of the parent compound; see the reference for this entry for the naming 
    convention applied to compounds other than pigments }
    \item{\code{polarity}:}{Object of class \code{"factor"}, ionization mode of 
    data in the database }
    \item{\code{num_entries}:}{Object of class \code{"integer"}, number of total
    entries (adducts) in the database }
    \item{\code{num_compounds}:}{Object of class \code{"integer"}, number of 
    parent compounds represented in the database (should be < 
    \code{num_entries})}
  }
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
  \code{\linkS4class{LOBdbase-class}},
  \code{\link{LOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{doLOBscreen}},
  \code{\link{default.LOBdbase}},
  \code{\link{default.componentCompTable}},
  \code{\link{default.adductHierarchies}},
  \code{\link{default.acylRanges}},
  \code{\link{default.oxyRanges}}
}
\examples{
## generate the default positive ionization mode database

\donttest{
LOBdbase.pos = generateLOBdbase(polarity = "positive", gen.csv = FALSE,
                                component.defs = NULL, AIH.defs = NULL,
                                acyl.ranges = NULL, oxy.ranges = NULL)
}
}
\keyword{methods}
