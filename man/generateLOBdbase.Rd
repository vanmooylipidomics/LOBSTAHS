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

For each lipid class or compound specified in the component definitions table,
the field \code{DB_gen_compound_type} must contain one of five values:
"DB_acyl_iteration," "DB_unique_species," "basic_component," "adduct_pos,"
or "adduct_neg." The last three compound types are reserved for definition of
basic components such as acteonitrile or acetate and for definition of adduct
ion types; new entries of these types should only be created in the  compound
table when a new adduct or basic component must be specified. The first two
compound types are used to define the way \code{generateLOBdbase} creates its
databases. There are essentially two ways \code{generateLOBdbase} creates
database entries in LOBSTAHS.

For compounds of \code{DB_gen_compound_type} = "DB_unique_species" (the simpler
of the two cases), database entries will be created only for adduct ions of the
single compound specified. The latter type should be used for pigments and other
lipids that do not have acyl groups, or when the user does not wish to consider
any possible variation in acyl properties. In this case, the exact mass of the
complete (neutral) molecule should be specified in the component definitions
(i.e., component composition) table.

Alternatively, for compounds of \code{DB_gen_compound_type} = 
"DB_acyl_iteration", \code{generateLOBdbase} will create database entries for
adduct ions of multiple molecular species within the lipid class based on the
ranges of acyl properties and oxidation states given for the class in 
\code{acyl.ranges} and \code{oxy.ranges} (see below). In this case, the compound
table should be used to define the exact mass of a "base fragment" for the lipid
class. Using this "base fragment" as a starting point, \code{generateLOBdbase}
creates multiple entries for molecules in the lipid class by iterative addition
of various combinations of fatty acids. In the case of IP-DAG and IP-MAG, the
base fragment includes the entire polar headgroup, the glycerol backbone, and
both carboxylic oxygen atoms in the fatty acid(s). In the case of TAG, the base
fragment is defined as the glycerol backbone plus the carboxylic oxygen atoms on
each of the three fatty acids. The base fragments for any new lipid classes for
which the user desires evaluation of a range of acyl properties should be
similarly defined.

Note that regardless of the \code{DB_gen_compound_type}, an adduct hierarchy
\strong{must} be specified in the adduct ion hierarchy matrix (see below) for
each compound or compound class specified in the
\code{Adduct_hierarchy_lookup_class} field of the component definitions table.
}
  \item{AIH.defs}{
File path to a .csv file containing empirical adduct ion hierarchy (AIH) data
for various pigments, lipids, and lipid classes. If nothing is specified,
\code{generateLOBdbase} will use the default AIH data (\code{\link{default.adductHierarchies}}).
Each compound or compound class for which there is an entry in the AIH
definitions table should have at least one corresponding entry in the
\code{Adduct_hierarchy_lookup_class} field of the component definitions table
(default, \code{\link{default.componentCompTable}}).

A Microsoft Excel spreadsheet template included with the package at
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_adduct_ion_hierarchies.xlsx can be
used to generate a .csv file with additional (or alternative) adduct hierarchy
data in a format appropriate for \code{generateLOBdbase}. Alternatively, the
spreadsheet may be downloaded at \url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_adduct_ion_hierarchies.xlsx}.
}
  \item{acyl.ranges}{
File path to a .csv file containing ranges of values for the total number of
acyl (i.e., fatty acid) carbon atoms to be considered during the
\emph{in silico} simulation of any lipid classses for which
\code{DB_gen_compound_type} is specified as "DB_acyl_iteration" in the component
definitions table, above. These include intact polar diacylglycerols (IP-DAG),
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
additional oxygen atoms to be considered during the \emph{in silico} simulation
of any lipid classses for which \code{DB_gen_compound_type} is specified as
"DB_acyl_iteration" in the component definitions table, above. If nothing is
specified, \code{generateLOBdbase} will use the default oxidation state ranges
in (\code{\link{default.oxyRanges}}).

A Microsoft Excel spreadsheet template included with the package at
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_addl_oxy_ranges.xlsx can be used to
generate a .csv file with custom oxidation state ranges in a format appropriate
for \code{generateLOBdbase}. Alternatively, the spreadsheet may be downloaded at
\url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_addl_oxy_ranges.xlsx}.
By default, \code{generateLOBdbase} considers 0-4 additional oxygen atoms on
each chemically possible IP-DAG, TAG, PUA, and FFA.
}
}
\details{Using the default structural property inputs described here,
\code{generateLOBdbase} can produce databases with entries for a wide range of
intact polar diacylglycerols (IP-DAG), triacylglycerols (TAG), polyunsaturated
aldehydes (PUAs), free fatty acids (FFA), and common photosynthetic pigments.
The default databases (as of January 2017) contain data on 18,067 and 15,404
unique compounds that can be identifed in positive and negative ion mode
spectra, respectively.

\strong{Note} that the default databases have been pre-generated (see
\code{\link{default.LOBdbase}}) and it is therefore unnecessary to call
\code{generateLOBdbase} with the default parameters.
}
\value{
  A \code{"\link{LOBdbase-class}"} object with the structure:
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
of oxidative stress biomarkers. \emph{Analytical Chemistry} 88:7154-7162 

Data for lipid classes BLL, PDPT, vGSL, sGSL, hGSL, hapGSL, and hapCER are as
described in Hunter J. E., M. J. Frada, H. F. Fredricks, A. Vardi, and B. A. S.
Van Mooy. 2015. Targeted and untargeted lipidomics of *Emiliania huxleyi* viral
infection and life cycle phases highlights molecular biomarkers of infection, 
susceptibility, and ploidy. \emph{Frontiers in Marine Science} 2:81
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\linkS4class{LOBdbase}},
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
