\name{LOBSet-class}
\Rdversion{1.1}
\docType{class}
\alias{LOBSet-class}
\alias{LOBSet}
\alias{show,LOBSet-method}

\title{Class \code{"LOBSet"}, Peak Data with Annotations, Isomers, and Compound Assignments}
\description{
A class for HPLC-MS peak data that have been screened and annotated using \code{LOBSTAHS} function \code{\link{doLOBscreen}}.
}
\section{Objects from the Class}{
A call of the form \code{new("LOBSet", ...)} will create a new, empty \code{"LOBdbase"} object with the slots indicated below. 
}
\section{Slots}{
  \describe{
    \item{\code{peakdata}:}{Object of class \code{"data.frame"}, containing peakdata by compound assignment }
    \item{\code{iso.C3r}:}{Object of class \code{"list"}, a list of the \code{match_ID}s of possible regioisomers of each compound }
    \item{\code{iso.C3f}:}{Object of class \code{"list"}, a list of the \code{match_ID}s of possible functional structural isomers of each compound }
    \item{\code{iso.C3c}:}{Object of class \code{"list"}, a list of the \code{match_ID}s of possible isobars of this compound }
    \item{\code{LOBscreen.diagnostics}:}{Object of class \code{"data.frame"}, numbers of peaks, peakgroups, adducts, and unique parent compounds present in the dataset after application of each screening criterion in \code{\link{doLOBscreen}} }
    \item{\code{LOBisoID.diagnostics}:}{Object of class \code{"data.frame"}, numbers of peakgroups and parent compounds to which the various isomer annotations have been assigned }
    \item{\code{LOBscreen.settings}:}{Object of class \code{"list"}, captures the settings used in \code{\link{doLOBscreen}} to generate the \code{"LOBSet"} }
    \item{\code{polarity}:}{Object of class \code{"factor"}, polarity of data in the \code{"LOBSet"} }
    \item{\code{sampnames}:}{Object of class \code{"character"}, number of samples represented in the \code{"LOBSet"} data }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "LOBSet")}: ... }
	 }
}
\references{
\item{Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, \dQuote{Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets}
}
\author{
James Collins, \email{james.r.collins@aya.yale.edu}
}
\seealso{
  \code{\link{doLOBscreen}},
  \code{\link{getLOBpeaklist}},
  \code{\linkS4class{xsAnnotate}}
}
\examples{
## return object information

showClass("LOBSet")
}
\keyword{classes}