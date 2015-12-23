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
\code{LOBSTAHS} Default Databases and Database Generation Parameters 
}
\description{
Default databases (\code{default.LOBdbase}) and .RData files containing the default ranges of structural properties used by \code{\link{generateLOBdbase}} to generate these databases. Also described is \code{default.rt.windows}, which contains the default retention time windows by lipid class used by \code{\link{doLOBscreen}} when \code{rt.restrict = TRUE}.
}
\usage{data(default.LOBdbase)}
\usage{data(default.adductHierarchies)}
\usage{data(default.acylRanges)}
\usage{data(default.oxyRanges)}
\usage{data(default.componentCompTable)}
\usage{data(default.rt.windows)}
\format{
  \describe{
    \item{\code{default.LOBdbase}}{A list of two \code{"\link{LOBdbase}"} objects, which are the default LOBSTAHS databases for positive and negative ion mode species, respectively. These were generated using \code{\link{generateLOBdbase}} with the default values defined in \code{default.adductHierarchies}, \code{default.acylRanges}, \code{default.oxyRanges}, and \code{default.componentCompTable}.}
    \item{\code{default.adductHierarchies}}{A data frame containing empirically-determined adduct ion hierarchy data, by lipid class.}
    \item{\code{default.acylRanges}}{A data frame containing ranges of numbers of acyl carbon atoms for which \strong{in silico} data are generated for each lipid class by \code{\link{generateLOBdbase}}.}
    \item{\code{default.oxyRanges}}{A data frame containing ranges of additional oxygen atoms to be considered on species of each lipid class when databases are generated with \code{\link{generateLOBdbase}}.}
    \item{\code{default.componentCompTable}}{A data frame that defines the elemental compositions of the various adducts, parent lipid "backbones," and pigments that are used by \code{\link{generateLOBdbase}}.}
    \item{\code{default.rt.windows}}{A data frame containing the default retention time data for various lipids and parent lipid classes that are used by \code{\link{doLOBscreen}} when \code{rt.restrict = TRUE}. These retention time windows are specific to the HPLC-MS method currently used in the Van Mooy Lab at Woods Hole Oceanographic Institution, where LOBSTAHS was developed. As described in \code{\link{doLOBscreen}}, users outside the Van Mooy Lab should supply their own retention time data.}
  }
  }
\details{
Empirical determination of the retention time window data in \code{default.rt.windows} and adduct ion hierarchies in \code{default.adductHierarchies} is described in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}. The default ranges for the structural properties given in the other files were chosen to yield databases that encompass a broad variety of moieites across lipid types. Microsoft Excel spreadsheet templates are included with the package in \url{../data-raw/xlsx/} for users wishing to modify any of the default data inputs. These templates can be used to generate .csv files in formats appropriate for \code{\link{generateLOBdbase}} and \code{\link{doLOBscreen}}.
}
\source{\url{http://github.com/vanmooylipidomics/LOBSTAHS/}
}
\references{
\item{Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
}
\seealso{
  \code{\link{doLOBscreen}},
  \code{\link{generateLOBdbase}},
  \code{\linkS4class{LOBdbase}}
}
\keyword{datasets}
