\name{doLOBscreen}
\alias{doLOBscreen}
\title{
Screen, annotate, and identify compounds in an xsAnnotate object
}
\description{
Primary function of the \code{LOBSTAHS} package. Screen, annotate, and assign
compound identities to peak groups in a CAMERA
\code{\link[CAMERA]{xsAnnotate}} object containing HPLC-ESI-MS lipid data.
Identify and annotate possible isomers and isobars.
}
\usage{
doLOBscreen(xsA, polarity = NULL, database = NULL, remove.iso = TRUE,
            rt.restrict =  TRUE, rt.windows = NULL, exclude.oddFA = TRUE,
            exclude.oxyrule =TRUE, match.ppm = NULL, retain.unidentified = TRUE)
}
\arguments{
  \item{xsA}{
A CAMERA \code{\link{xsAnnotate-class}} object with identified pseudospectra.
Must contain isotope data from \code{\link[CAMERA]{findIsotopes}} if user
elects \code{remove.iso = TRUE}).
}
  \item{polarity}{
Specify polarity mode of data in \code{xsA}: "positive" or "negative."
\code{doLOBscreen} will attempt to detect the polarity mode if nothing is given.
}
  \item{database}{
Specify the \code{LOBdbase} from which compound identities and adduct ion
hierarchy data are to be applied. If nothing is specified, \code{doLOBscreen}
will use the LOBSTAHS default database (\code{\link{default.LOBdbase}}) for the
appropriate ionization mode.

The database generation function \code{\link{generateLOBdbase}} can be used to
create a LOBdbase from structural property ranges specified by the user in a
series of input tables. \code{\link{loadLOBdbase}} can be used to import and
reconstruct a \code{"LOBdbase"} object from a previously generated database that
was saved as a .csv file.
}
  \item{remove.iso}{
Should secondary isotope peaks be removed? (If \code{TRUE}, \code{xsA} must
contain isotope data obtained using \code{\link[CAMERA]{findIsotopes}}.)
}
  \item{rt.restrict}{
Should lipid class retention time screening criteria be applied to putative
compound assignments?
}
  \item{rt.windows}{
File path to a .csv file containing retention time "window" data to be used for
screening by lipid class, if \code{rt.restrict = TRUE}. If nothing is specified,
\code{doLOBscreen} will use the package default windows
(\code{\link{default.rt.windows}}), which are specific to the HPLC-MS method
currently used in the Van Mooy Lab at Woods Hole Oceanographic Institution.

Because the user is advised to provide retention time data specific to his/her
HPLC-MS method, failure to specify a value for \code{rt.windows} will result in
a warning. A Microsoft Excel spreadsheet template included with the package at
Resources/library/LOBSTAHS/doc/xlsx/LOBSTAHS_rt.windows.xlsx can be
used to generate a .csv file of retention time data in a format appropriate for
\code{rt.windows}. Alternatively, the spreadsheet file may be downloaded
\href{https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/inst/doc/xlsx/LOBSTAHS_rt.windows.xlsx}{from the package GitHub repository}.
}
  \item{exclude.oddFA}{
Should compound assignments with an odd total number of acyl carbon atoms be
excluded? (Applies only to assignments where the parent lipid class is TAG,
IP-DAG, PUA, FFA or IP-MAG.) Useful if data are (or are believed to be) of exclusively
eukaryotic origin.
}
  \item{match.ppm}{
\emph{m/z} tolerance (in ppm) used for matching observed data against the
calculated \emph{m/z}'s in the database
}
  \item{exclude.oxyrule}{
Should compound assignments include physically unlikey fatty acid molecules where the total number of double
bonds are less than the oxidation state of the molecule? (Applies only to assignments where the parent lipid 
class is TAG, IP-DAG, PUA, FFA or IP-MAG.) Useful if such rare molecules are likely to be found in sample.
}
  \item{retain.unidentified}{
Should data for unidentified or discarded features be retained and appended to 
user's results? Useful if untargeted follow-on data analysis is anticipated, or
when user wants to retain data for \emph{all} features present in the original
dataset, not just those for which a LOBSTAHS identity was found.
}
}
\details{\code{doLOBscreen} draws compound identities from a \code{LOBdbase}
database. The function applies various retention time and adduct ion hierarchy
screening criteria to winnow the list of putative compound assignments. It
returns a table of annotated peak data with compound assignments and various
annotation codes to assist the user in interpretation and follow-on data
analysis.
}
\value{
  A \code{"\link{LOBSet-class}"} object.
}
\references{
The LOBSTAHS package is presented in:

Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy. 2016.
LOBSTAHS: An adduct-based lipidomics strategy for discovery and identification
of oxidative stress biomarkers. \emph{Analytical Chemistry} 88:7154-7162,
doi:\href{http://dx.doi.org/10.1021/acs.analchem.6b01260}{10.1021/acs.analchem.6b01260}.

Data for lipid classes BLL, PDPT, vGSL, sGSL, hGSL, hapGSL, and hapCER are as
described in:

Hunter J. E., M. J. Frada, H. F. Fredricks, A. Vardi, and B. A. S.
Van Mooy. 2015. Targeted and untargeted lipidomics of \emph{Emiliania huxleyi}
viral infection and life cycle phases highlights molecular biomarkers of
infection, susceptibility, and ploidy. \emph{Frontiers in Marine Science} 2:81,
doi:\href{http://dx.doi.org/10.3389/fmars.2015.00081}{10.3389/fmars.2015.00081}

Fulton, J. M., H. F. Fredricks, K. D. Bidle, A. Vardi, B. J. Kendrick, G. R.
DiTullio, and B. A. S. Van Mooy. 2014. Novel molecular determinants of viral
susceptibility and resistance in the lipidome of \emph{Emiliania huxleyi},
\emph{Environmental Microbiology} 16(4):1137-1149,
doi:\href{http://dx.doi.org/10.1111/1462-2920.12358}{10.1111/1462-2920.12358}.
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
## screen & annotate xsAnnotate object from the PtH2O2lipids dataset using all
## screening options, while retaining unidentified/discarded features
library(PtH2O2lipids)

\donttest{
myPtH202LOBSet = doLOBscreen(ptH2O2lipids$xsAnnotate, polarity = "positive",
                             database = NULL, remove.iso = TRUE,
                             rt.restrict =  TRUE, rt.windows = NULL,
                             exclude.oddFA = TRUE, match.ppm = 2.5, exclude.oxyrule =TRUE,
                             retain.unidentified = TRUE)

## show some diagnostics
LOBscreen_diagnostics(myPtH202LOBSet) # screening diagnostics
LOBisoID_diagnostics(myPtH202LOBSet) # isomer identification data
}
}
\keyword{methods}
