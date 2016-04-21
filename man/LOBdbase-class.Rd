\name{LOBdbase-class}
\Rdversion{1.1}
\docType{class}
\alias{LOBdbase-class}
\alias{LOBdbase}
\alias{[,LOBdbase,ANY,ANY,ANY-method}
\alias{show,LOBdbase-method}

\title{Class LOBdbase: A class for LOBSTAHS lipid-oxylipin databases}
\description{
A class for \code{LOBSTAHS} databases that contain a combination of \emph{in silico} and empirical data for a wide range of lipids, oxidized lipids, and oxylipins. A \code{"LOBdbase"} object can be generated using \code{\link{generateLOBdbase}}, or imported from a .csv file of proper format using \code{\link{loadLOBdbase}}.
}
\section{Objects from the Class}{
A call of the form \code{new("LOBdbase", ...)} will create a new, empty \code{"LOBdbase"} object with the slots indicated below. 
}
\section{Slots}{
  \describe{
    \item{\code{frag_ID}:}{Object of class \code{"integer"}, a unique identifier for this molecular species }
    \item{\code{mz}:}{Object of class \code{"numeric"}, the calculated \emph{m/z} of this species }
    \item{\code{exact_parent_neutral_mass}:}{Object of class \code{"numeric"}, the calculated (monoisotopic) exact mass of the parent compound of this species }
    \item{\code{lipid_class}:}{Object of class \code{"factor"}, the parent lipid class of this species }
    \item{\code{species}:}{Object of class \code{"character"}, the lipid subclass }
    \item{\code{adduct}:}{Object of class \code{"factor"}, the adduct ion represented by this entry }
    \item{\code{adduct_rank}:}{Object of class \code{"integer"}, the relative abundance ranking of this adduct relative to the other adducts of the same parent compound }
    \item{\code{FA_total_no_C}:}{Object of class \code{"integer"}, total number of acyl (fatty acid) carbon atoms in the parent compound; \code{NA} if \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{FA_total_no_DB}:}{Object of class \code{"integer"}, total number of acyl (fatty acid) carbon-carbon double bonds in the parent compound; \code{NA} if \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{degree_oxidation}:}{Object of class \code{"integer"}, number of additional oxygen atoms present }
    \item{\code{parent_elem_formula}:}{Object of class \code{"character"}, elemental formula of the parent compound }
    \item{\code{parent_compound_name}:}{Object of class \code{"character"}, name of the parent compound; see the reference for this entry for the naming convention applied to compounds other than pigments }
    \item{\code{polarity}:}{Object of class \code{"factor"}, ionization mode of data in the database }
    \item{\code{num.entries}:}{Object of class \code{"integer"}, number of total entries (adducts) in the database }
    \item{\code{num.compounds}:}{Object of class \code{"integer"}, number of parent compounds represented in the database (should be < \code{num.entries})}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "LOBdbase")}: ... }
    \item{show}{\code{signature(object = "LOBdbase")}: ... }
	 }
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\link{generateLOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{doLOBscreen}}
}
\examples{
## return object information

showClass("LOBdbase")
}
\keyword{classes}
