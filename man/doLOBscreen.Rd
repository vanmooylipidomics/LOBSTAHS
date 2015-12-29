\name{doLOBscreen}
\alias{doLOBscreen}
\title{
Screen, Annotate, and Identify Compounds in an xsAnnotate Object
}
\description{
Primary function of the \code{LOBSTAHS} package. Screen, annotate, and assign compound identities to peak groups in a CAMERA \code{"\link[CAMERA]{xsAnnotate}"} object containing HPLC-ESI-MS lipid data. Identify and annotate possible isomers and isobars.
}
\usage{
doLOBscreen(xsA, polarity = NULL, database = NULL, remove.iso = TRUE,
            rt.restrict =  TRUE, rt.windows = NULL, exclude.oddFA = TRUE,
            match.ppm = NULL)
}
\arguments{
  \item{xsA}{
A CAMERA \code{"xsAnnotate"} object with identified pseudospectra. Must contain isotope data from \code{\link[CAMERA]{findIsotopes}} if user elects \code{remove.iso = TRUE}).
}
  \item{polarity}{
Specify polarity mode of data in \code{xsA}: "positive" or "negative." \code{doLOBscreen} will attempt to detect the polarity mode if nothing is given.
}
  \item{database}{
Specify the \code{LOBdbase} from which compound identities and adduct ion hierarchy data are to be applied. If nothing is specified, \code{doLOBscreen} will use the LOBSTAHS default database (\code{\link{default.LOBdbase}}) for the appropriate ionization mode.

The database generation function \code{\link{generateLOBdbase}} can be used to create a LOBdbase from structural property ranges specified by the user in a series of input tables. \code{\link{loadLOBdbase}} can be used to import and reconstruct a \code{"LOBdbase"} object from a previously generated database that was saved as a .csv file.
}
  \item{remove.iso}{
Should secondary isotope peaks be removed? (If \code{TRUE}, \code{xsA} must contain isotope data obtained using \code{\link[CAMERA]{findIsotopes}}.)
}
  \item{rt.restrict}{
Should lipid class retention time screening criteria be applied to putative compound assignments?
}
  \item{rt.windows}{
File path to a .csv file containing retention time "window" data to be used for screening by lipid class, if \code{rt.restrict = TRUE}. If nothing is specified, \code{doLOBscreen} will use the package default windows (\code{\link{default.rt.windows}}), which are specific to the HPLC-MS method currently used in the Van Mooy Lab at Woods Hole Oceanographic Institution.

Because the user is advised to provide retention time data specific to his/her HPLC-MS method, failure to specify a value for \code{rt.windows} will result in a warning. A Microsoft Excel spreadsheet template included with the package at Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_lipid_class_rt_windows.xlsx can be used to generate a .csv file of retention time data in a format appropriate for \code{rt.windows}. Alternatively, the spreadsheet file may be downloaded at \url{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_lipid_class_rt_windows.xlsx}.
}
  \item{exclude.oddFA}{
Should compound assignments with an odd total number of acyl carbon atoms be excluded? (Applies only to assignments where the parent lipid class is TAG, IP-DAG, PUA, or FFA.) Useful if data are (or are believed to be) of exclusively eukaryotic origin.
}
  \item{match.ppm}{
\emph{m/z} tolerance (in ppm) used for matching observed data against the calculated \emph{m/z}'s in the database
}
}
\details{\code{doLOBscreen} draws compound identities from a \code{LOBdbase} database. The function applies various retention time and adduct ion hierarchy screening criteria to winnow the list of putative compound assignments. It returns a table of annotated peak data with compound assignments and various annotation codes to assist the user in interpretation and follow-on data analysis.
}
\value{
  A \code{"\link{LOBSet}"} object.
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\linkS4class{LOBSet}},
  \code{\linkS4class{LOBdbase}},
  \code{\link{loadLOBdbase}},
  \code{\link{getLOBpeaklist}},
  \code{\link{generateLOBdbase}},
  \code{\link{default.LOBdbase}},
  \code{\link{default.rt.windows}},
  \code{\linkS4class{xcmsSet}},
  \code{\linkS4class{xsAnnotate}},
  \code{\link[CAMERA]{findIsotopes}}
}
\examples{
## screen & annotate "xsAnnotate" object from the PtH2O2lipids dataset using all
## screening options
library(LOBSTAHS)
library(PtH2O2lipids)

\donttest{
myPtH202LOBSet = doLOBscreen(ptH2O2lipids$xsAnnotate, polarity = "positive",
                             database = NULL, remove.iso = TRUE,
                             rt.restrict =  TRUE, rt.windows = NULL,
                             exclude.oddFA = TRUE, match.ppm = 2.5)

## show some diagnostics
myPtH202LOBSet@LOBscreen.diagnostics # screening diagnostics
myPtH202LOBSet@LOBisoID.diagnostics # isomer identification data
}
}
\keyword{methods}
