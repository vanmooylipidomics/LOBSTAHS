\name{LOBdefaults}
\alias{LOBdefaults}
\alias{default.LOBdbase}
\alias{default.adductHierarchies}
\alias{default.acylRanges}
\alias{default.oxyRanges}
\alias{default.componentCompTable}
\alias{default.rt.windows}
\docType{data}
\title{
LOBSTAHS default databases and database generation parameters 
}
\description{
Default databases (in \code{default.LOBdbase}), and the .RData files containing 
default ranges of structural properties used by \code{\link{generateLOBdbase}} 
to generate these databases.

Also described is \code{default.rt.windows}, which contains the default 
retention time windows (by lipid class) used by \code{\link{doLOBscreen}} when 
\code{rt.restrict = TRUE}.
}
\usage{data(default.LOBdbase)
data(default.adductHierarchies)
data(default.acylRanges)
data(default.oxyRanges)
data(default.componentCompTable)
data(default.rt.windows)}
\format{
  \describe{
    \item{\code{default.LOBdbase}}{A list of two \code{"\link{LOBdbase}"} 
    objects, which are the default LOBSTAHS databases for positive and negative 
    ion mode species, respectively. These were generated using 
    \code{\link{generateLOBdbase}} with the default values defined in 
    \code{default.adductHierarchies}, \code{default.acylRanges}, 
    \code{default.oxyRanges}, and \code{default.componentCompTable}.}
    \item{\code{default.adductHierarchies}}{A data frame containing 
    empirically-determined adduct ion hierarchy data, by lipid class.}
    \item{\code{default.acylRanges}}{A data frame containing ranges of numbers 
    of acyl carbon atoms for which \emph{in silico} data are generated for each 
    lipid class by \code{\link{generateLOBdbase}}.}
    \item{\code{default.oxyRanges}}{A data frame containing ranges of additional
    oxygen atoms to be considered on species of each lipid class when databases 
    are generated with \code{\link{generateLOBdbase}}.}
    \item{\code{default.componentCompTable}}{A data frame that defines the 
    elemental compositions of the various adducts, parent lipid "backbones," and
    pigments that are used by \code{\link{generateLOBdbase}}.}
    \item{\code{default.rt.windows}}{A data frame containing the default 
    retention time data for various lipids and parent lipid classes that are 
    used by \code{\link{doLOBscreen}} when \code{rt.restrict = TRUE}. These 
    retention time windows are specific to the HPLC-MS method currently used in 
    the Van Mooy Lab at Woods Hole Oceanographic Institution, where LOBSTAHS was
    developed. As described in \code{\link{doLOBscreen}}, users outside the Van 
    Mooy Lab should supply their own retention time data.}
  }
  }
\details{
Empirical determination of the retention time window data in 
\code{default.rt.windows} and adduct ion hierarchies in 
\code{default.adductHierarchies} is described in the reference below.

The default ranges for the structural properties given in the other files were 
chosen to yield databases that encompass a broad variety of moieites across 
lipid types. Microsoft Excel spreadsheet templates are included with the package
in Resources/library/LOBSTAHS/doc for users wishing to modify any of the default
data inputs. Alternatively, the spreadsheet files may be downloaded from 
\url{https://github.com/vanmooylipidomics/LOBSTAHS/tree/master/inst/doc/xlsx}. 
These templates can be used to generate .csv files in formats appropriate for 
\code{\link{generateLOBdbase}} and \code{\link{doLOBscreen}}.
}
\source{\url{http://github.com/vanmooylipidomics/LOBSTAHS/}
}
\value{
  Various list and data.frame objects (as indicated above).
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy. 2016.
LOBSTAHS: An adduct-based lipidomics strategy for discovery and identification
of oxidative stress biomarkers. Analytical Chemistry.
}
\seealso{
  \code{\link{doLOBscreen}},
  \code{\link{generateLOBdbase}},
  \code{\linkS4class{LOBdbase}}
  \code{\link{LOBdbase}}
}
\keyword{datasets}
