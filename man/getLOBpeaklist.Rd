\name{getLOBpeaklist}
\alias{getLOBpeaklist}
\title{
Export Screened \code{LOBSTAHS} Peaklist with Compound Assignments
}
\description{
Extracts screened peak data, compound assignments, annotation codes, and (optionally) isomer data from a \code{\link[LOBSTAHS]{"LOBSet"}} object. Returns a table containing the requested data, with option to export to a .csv file.
}
\usage{
getLOBpeaklist(LOBSet, include.iso = TRUE, gen.csv = FALSE)
}
\arguments{
  \item{LOBSet}{
A LOBSTAHS \code{"LOBSet"} object.
}
  \item{include.iso}{
Specify whether isomer and isobar identification data should be included in peaklist.
}
  \item{gen.csv}{
Should a .csv file be generated in addition to the data frame that is returned?
}
}
\details{\code{getLOBpeaklist} extracts data from all slots in a given \code{"LOBSet"} object and flows it into a data frame. Annotation codes indicating compliance with the adduct ion hierarchy screening criteria are appended to the \strong{m/z}, retention time, and peak area data for each xcms peakgroup remaining in the final dataset. If \code{include.iso = TRUE}, three additional columns containing the \code{match_IDs} of the possible isomers and isobars for each peakgroup are also appended. \code{getLOBpeaklist} does not export any of the diagnostic data in the \code{LOBSet} \code{@LOBscreen.diagnostics} or \code{@LOBisoID.diagnostics} slots.
}
\value{
  A data frame with the following structure:
  \samp{
  \code{match_ID}           integer; a unique identifier for each assignment of
                            a compound to a peakgroup (multiple match_IDs can
                            exist for a peakgroup if the group was assigned
                            multiple possible compound identities
  \code{compound_name}      character; name of the compound; see reference for
                            naming convention applied to compounds other than
                            pigments
  \code{elem_formula}       character; elemental formula of the compound
  \code{LOBdbase.mz}        numeric; calculated monoisotopic \strong{m/z} of the
                            compound from database
  \code{peakgroup.mz}       numeric; mean observed \strong{m/z} of the feature
                            across all samples in which it was identified
  \code{LOBdbase.ppm.match} numeric; ppm deviation between observed and
                            theoretical \strong{m/z}
  \code{peakgroup.rt}       numeric; mean observed retention time of the feature
                            across all samples in which it was identified
  \code{peakgroup.mzmin}    numeric; minimum observed \strong{m/z} across
                            samples
  \code{peakgroup.mzmax}    numeric; maximum observed \strong{m/z} across
                            samples
  \code{peakgroup.rtmin}    numeric; minimum observed retention time across
                            samples
  \code{peakgroup.rtmax}    numeric; maximum observed retention time across
                            samples

  several columns (numeric) containing integrated peak area data for this group
  (one column for each sample in the dataset)

  \code{xcms_peakgroup}     integer; the xcms peakgroup identifier
  \code{CAMERA_pseudospectrum}
                            integer; the CAMERA pseudospectrum identifier
  \code{LOBdbase.frag_ID}   integer; the LOBdbase fragment ID corresponding to
                            the adduct of this compound for which data are 
                            reported; this is the dominant adduct of the
                            compound, according to the adduct ion hierarchy
                            rules for the parent lipid class
  \code{LOBdbase.exact_parent_neutral_mass}
                            numeric; the calculated exact (monoisotopic) mass of
                            the compound
  \code{lipid_class}        character; parent lipid class of this compound
  \code{species}            character; if a pigment or IP-DAG, the specific
                            compound species
  \code{major.adduct}       character; adduct of the compound for which data is
                            reported
  \code{FA_total_no_C}      integer; total number of acyl (fatty acid) carbon
                            atoms in this compound; \code{NA} if \code{lipid_class}
                            is not TAG, IP-DAG, PUA, or FFA
  \code{FA_total_no_DB}     integer; total number of acyl (fatty acid)
                            carbon-carbon double bonds in this compound; \code{NA}
                            if \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA
  \code{degree_oxidation}   integer; number of additional oxygen atoms present
  \code{C1}-\code{C6b}      columns (integer) containing binary indicators for
                            each possible annotation code applied by \code{doLOBscreen}
  \code{casecodes}          character string containing list of all codes
                            applied to this assignment
                            
  and, if \code{include.iso = TRUE}:
  
  \code{iso.C3r.match_ID}   character string of integers; the \code{match_ID}s of all
                            possible regioisomers of this compound
  \code{iso.C3f.match_ID}   character string of integers; the \code{match_ID}s of all
                            possible functional structural isomers of this
                            compound
  \code{iso.C3c.match_ID}   character string of integers; the \code{match_ID}s of all
                            possible isobars of this compound (i.e., those
                            assignments which cannot be distinguished from this
                            one given the applied ppm tolerance)
  }
}
\references{
\item{Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\linkS4class[LOBSTAHS]{LOBSet}},
  \code{\link[LOBSTAHS]{doLOBscreen}}
}
\examples{
## export peaklist, with isomer data
library(xcms)
library(CAMERA)
library(LOBSTAHS)
library(PtH2O2lipids)

PtH2O2.peakdata = getLOBpeaklist(PtH2O2lipids$LOBSet, include.iso = TRUE,
                                 gen.csv = FALSE)
}
\keyword{methods}
