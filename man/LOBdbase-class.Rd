\name{LOBdbase-class}
\Rdversion{1.1}
\docType{class}
\alias{LOBdbase-class}
\alias{show,LOBdbase-method}
\alias{frag_ID,LOBdbase-method}
\alias{frag_ID}
\alias{frag_ID<-,LOBdbase-method}
\alias{frag_ID<-}
\alias{mz,LOBdbase-method}
\alias{mz}
\alias{mz<-,LOBdbase-method}
\alias{mz<-}
\alias{exact_parent_neutral_mass,LOBdbase-method}
\alias{exact_parent_neutral_mass}
\alias{exact_parent_neutral_mass<-,LOBdbase-method}
\alias{exact_parent_neutral_mass<-}
\alias{lipid_class,LOBdbase-method}
\alias{lipid_class}
\alias{lipid_class<-,LOBdbase-method}
\alias{lipid_class<-}
\alias{species,LOBdbase-method}
\alias{species}
\alias{species<-,LOBdbase-method}
\alias{species<-}
\alias{adduct,LOBdbase-method}
\alias{adduct}
\alias{adduct<-,LOBdbase-method}
\alias{adduct<-}
\alias{adduct_rank,LOBdbase-method}
\alias{adduct_rank}
\alias{adduct_rank<-,LOBdbase-method}
\alias{adduct_rank<-}
\alias{FA_total_no_C,LOBdbase-method}
\alias{FA_total_no_C}
\alias{FA_total_no_C<-,LOBdbase-method}
\alias{FA_total_no_C<-}
\alias{FA_total_no_DB,LOBdbase-method}
\alias{FA_total_no_DB}
\alias{FA_total_no_DB<-,LOBdbase-method}
\alias{FA_total_no_DB<-}
\alias{degree_oxidation,LOBdbase-method}
\alias{degree_oxidation}
\alias{degree_oxidation<-,LOBdbase-method}
\alias{degree_oxidation<-}
\alias{parent_elem_formula,LOBdbase-method}
\alias{parent_elem_formula}
\alias{parent_elem_formula<-,LOBdbase-method}
\alias{parent_elem_formula<-}
\alias{parent_compound_name,LOBdbase-method}
\alias{parent_compound_name}
\alias{parent_compound_name<-,LOBdbase-method}
\alias{parent_compound_name<-}
\alias{polarity,LOBdbase-method}
\alias{polarity<-,LOBdbase-method}
\alias{num_entries,LOBdbase-method}
\alias{num_entries}
\alias{num_entries<-,LOBdbase-method}
\alias{num_entries<-}
\alias{num_compounds,LOBdbase-method}
\alias{num_compounds}
\alias{num_compounds<-,LOBdbase-method}
\alias{num_compounds<-}

\title{Class LOBdbase: A class for LOBSTAHS lipid-oxylipin databases}
\description{
A class for \code{LOBSTAHS} databases that contain a combination of 
\emph{in silico} and empirical data for a wide range of lipids, oxidized lipids,
and oxylipins.
}
\section{Objects from the class}{
Objects can be created using the simulation function 
\code{\link{generateLOBdbase}} (preferred; satisfies the needs of most
users), imported from a .csv file of proper format using 
\code{\link{loadLOBdbase}}, or created using the rudimentary constructor 
\code{\link{LOBdbase}} (least preferred). A 'LOBdbase' can also be created by 
calls of the form \code{new("LOBdbase", ...)}.
}
\section{Slots}{
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
\section{Methods}{
  \describe{
    \item{show}{
      \code{signature(object = "LOBdbase")}: ...
    }
    \item{polarity}{
      \code{signature(object = "LOBdbase")}: get \code{polarity} slot
    }
    \item{num_compounds}{
      \code{signature(object = "LOBdbase")}: get \code{num_compounds} slot
    }
    \item{num_entries}{
      \code{signature(object = "LOBdbase")}: get \code{num_entries} slot
    }
    \item{frag_ID}{
      \code{signature(object = "LOBdbase")}: get \code{frag_ID} slot
    }
    \item{exact_parent_neutral_mass}{
      \code{signature(object = "LOBdbase")}: get 
      \code{exact_parent_neutral_mass} slot
    }
    \item{lipid_class}{
      \code{signature(object = "LOBdbase")}: get \code{lipid_class} slot
    }
    \item{species}{
      \code{signature(object = "LOBdbase")}: get \code{species} slot
    }
    \item{adduct}{
      \code{signature(object = "LOBdbase")}: get \code{adduct} slot
    }
    \item{adduct_rank}{
      \code{signature(object = "LOBdbase")}: get \code{adduct_rank} slot
    }
    \item{FA_total_no_C}{
      \code{signature(object = "LOBdbase")}: get \code{FA_total_no_C} slot
    }
    \item{FA_total_no_DB}{
      \code{signature(object = "LOBdbase")}: get \code{FA_total_no_DB} slot
    }
    \item{degree_oxidation}{
      \code{signature(object = "LOBdbase")}: get \code{degree_oxidation} slot
    }
    \item{parent_elem_formula}{
      \code{signature(object = "LOBdbase")}: get \code{parent_elem_formula} slot
    }
    \item{parent_compound_name}{
      \code{signature(object = "LOBdbase")}: get \code{parent_compound_name} 
      slot
    }
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
  \code{\link{generateLOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{doLOBscreen}},
  \code{\link{LOBdbase}}
}
\examples{
## return object information

showClass("LOBdbase")
}
\keyword{classes}
