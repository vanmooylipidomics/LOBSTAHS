\name{getLOBpeaklist}
\alias{getLOBpeaklist}
\title{
Export screened LOBSTAHS peaklist with compound assignments
}
\description{
Extracts screened peak data, compound assignments, annotation codes, and 
(optionally) isomer data from a \code{\link{LOBSet-class}} object. Returns a 
table containing the requested data, with option to export to a .csv file.
}
\usage{
getLOBpeaklist(LOBSet, include.iso = TRUE, include.unidentified = TRUE,
               gen.csv = FALSE)
}
\arguments{
  \item{LOBSet}{
A LOBSTAHS \code{"LOBSet"} object.
}
  \item{include.iso}{
Specify whether isomer and isobar identification data should be included in 
peaklist.
}
  \item{include.unidentified}{
If data for unidentified or discarded features is present in the \code{LOBSet},
should these data be exported with the peaklist? A value of \code{FALSE} will
yield a peaklist that includes data only for those features to which one or more
LOBSTAHS compound identities have been assigned.}
  \item{gen.csv}{
Should a .csv file be generated in addition to the data frame that is returned?
}
}
\details{\code{getLOBpeaklist} extracts data from all slots in a given 
\code{"LOBSet"} object and flows it into a data frame. Annotation codes 
indicating compliance with the adduct ion hierarchy screening criteria are 
appended to the \emph{m/z}, retention time, and peak area data for each xcms 
peakgroup remaining in the final dataset.

If \code{include.iso = TRUE}, three additional columns containing the 
\code{match_IDs} of the possible isomers and isobars for each peakgroup are also
appended. \code{getLOBpeaklist} does not export any of the diagnostic data in 
the \code{LOBSet} \code{LOBscreen_diagnostics} or \code{LOBisoID_diagnostics} 
slots.

When the \code{LOBSet} includes peak data for features that were not identified
during screening or were discarded (i.e., if the user specified
\code{retain.unidentified = TRUE} when calling \code{\link{doLOBscreen}}), the
user is given the option to export these data when using
\code{getLOBpeaklist}. This is useful when untargeted follow-on data analysis is
anticipated, or the user simply wants to export data for \emph{all} features
present in the original dataset, not just those for which a LOBSTAHS identity
was found.
}
\value{
  A data frame with the following structure:
    \describe{
    \item{\code{match_ID}:}{Object of class \code{"integer"}, unique identifier 
    for each assignment of a compound to a peakgroup (multiple \code{match_ID}s 
    can exist for a peakgroup if the group was assigned multiple compound 
    identities }
    \item{\code{compound_name}:}{Object of class \code{"character"}, name of 
    compound; see reference for naming convention applied to compounds other 
    than pigments }
    \item{\code{elem_formula}:}{Object of class \code{"character"}, empirical 
    formula of compound }
    \item{\code{LOBdbase_mz}:}{Object of class \code{"numeric"}, calculated 
    \emph{m/z} of the adduct for which data in this group are reported; obtained
    database }
    \item{\code{peakgroup_mz}:}{Object of class \code{"numeric"}, mean observed
    \emph{m/z} of the feature in this peakgroup across all samples in which it 
    was identified }
    \item{\code{LOBdbase_ppm_match}:}{Object of class \code{"numeric"}, ppm 
    deviation between observed and calculated \emph{m/z} }
    \item{\code{peakgroup_rt}:}{Object of class \code{"numeric"}, mean observed
    retention time of the feature in this peakgroup across all samples in which 
    it was identified }
    \item{\code{peakgroup_mzmin}:}{Object of class \code{"numeric"}, minimum 
    observed \emph{m/z} of feature across samples }
    \item{\code{peakgroup_mzmax}:}{Object of class \code{"numeric"}, maximum 
    observed \emph{m/z} of feature across samples }
    \item{\code{peakgroup_rtmin}:}{Object of class \code{"numeric"}, minimum 
    observed retention time of feature across samples }
    \item{\code{peakgroup_rtmax}:}{Object of class \code{"numeric"}, maximum 
    observed retention time of feature across samples }
    \item{peak area data:}{Several objects of class \code{"numeric"}, containing
    integrated peak area data for this group by sample (one column for each 
    sample in the dataset) }
    \item{\code{xcms_peakgroup}:}{Object of class \code{"integer"}, the xcms 
    \code{xcmsSet} peakgroup identifier }
    \item{\code{CAMERA_pseudospectrum}:}{Object of class \code{"integer"}, the 
    CAMERA \code{xsAnnotate} pseudospectrum identifier }
    \item{\code{LOBdbase_frag_ID}:}{Object of class \code{"integer"}, the 
    \code{LOBdbase} fragment ID corresponding to the adduct of this compound for
    which data are  reported; this is the dominant adduct of the compound 
    according to the adduct ion hierarchy rules for the parent lipid class }
    \item{\code{LOBdbase_exact_parent_neutral_mass}:}{Object of class 
    \code{"numeric"}, the calculated exact (monoisotopic) mass of the compound; 
    from database }
    \item{\code{lipid_class}:}{Object of class \code{"factor"}, parent lipid 
    class of this compound }
    \item{\code{species}:}{Object of class \code{"character"}, if a pigment or 
    IP-DAG, the specific compound species }
    \item{\code{major_adduct}:}{Object of class \code{"factor"}, adduct of the 
    compound for which data in this entry is reported }
    \item{\code{FA_total_no_C}:}{Object of class \code{"integer"}, total number 
    of acyl (fatty acid) carbon atoms in this compound; \code{"NA"} if 
    \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{FA_total_no_DB}:}{Object of class \code{"integer"}, total number
    of acyl (fatty acid) carbon-carbon double bonds in this compound; 
    \code{"NA"} if \code{lipid_class} is not TAG, IP-DAG, PUA, or FFA }
    \item{\code{degree_oxidation}:}{Object of class \code{"integer"}, number of 
    additional oxygen atoms present on this compound, compared with its 
    unoxidized parent }
    \item{\code{C1}-\code{C6b}:}{Several objects of class \code{"integer"}, 
    containing binary indicators for each possible annotation code applied by 
    \code{doLOBscreen} }
    \item{\code{casecodes}:}{Object of class \code{"character"}, character 
    string containing list of all codes applied to this assignment }
    \item{\code{iso_C3r_match_ID}:}{Object of class \code{"character"}, 
    character string of integer containing the \code{match_ID}s of all possible 
    regioisomers of this compound (if \code{include.iso = TRUE}) }
    \item{\code{iso_C3f_match_ID}:}{Object of class \code{"character"}, 
    character string of integer containing the \code{match_ID}s of all possible 
    functional structural isomers of this compound (if \code{include.iso = TRUE}) }
    \item{\code{iso_C3c_match_ID}:}{Object of class \code{"character"}, 
    character string of integer containing the \code{match_ID}s of all possible 
    isobars of this compound (if \code{include.iso = TRUE}) }
  }
}
\references{
Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy. 2016.
LOBSTAHS: An adduct-based lipidomics strategy for discovery and identification
of oxidative stress biomarkers. \emph{Analytical Chemistry} 88:7154-7162 
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\link{LOBSet}},
  \code{\linkS4class{LOBSet}},
  \code{\link{doLOBscreen}}

}
\examples{
## export peaklist, with isomer data, but without data for unidentified features
library(PtH2O2lipids)

PtH2O2.peakdata = getLOBpeaklist(ptH2O2lipids$LOBSet, include.iso = TRUE,
                                 include.unidentified = FALSE, gen.csv = FALSE)
}
\keyword{methods}
