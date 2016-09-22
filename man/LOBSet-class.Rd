\name{LOBSet-class}
\Rdversion{1.1}
\docType{class}
\alias{LOBSet-class}
\alias{show,LOBSet-method}
\alias{LOBisoID_diagnostics,LOBSet-method}
\alias{LOBisoID_diagnostics}
\alias{LOBisoID_diagnostics<-,LOBSet-method}
\alias{LOBisoID_diagnostics<-}
\alias{peakdata,LOBSet-method}
\alias{peakdata}
\alias{peakdata<-,LOBSet-method}
\alias{peakdata<-}
\alias{polarity,LOBSet-method}
\alias{polarity}
\alias{polarity<-,LOBSet-method}
\alias{polarity<-}
\alias{sampnames,LOBSet-method}
\alias{sampnames}
\alias{sampnames<-,LOBSet-method}
\alias{sampnames<-}
\alias{LOBscreen_settings,LOBSet-method}
\alias{LOBscreen_settings}
\alias{LOBscreen_settings<-,LOBSet-method}
\alias{LOBscreen_settings<-}
\alias{LOBscreen_diagnostics,LOBSet-method}
\alias{LOBscreen_diagnostics}
\alias{LOBscreen_diagnostics<-,LOBSet-method}
\alias{LOBscreen_diagnostics<-}

\title{Class LOBSet: Peak data with annotations, isomers, and compound
assignments}
\description{
A class for HPLC-MS peak data that have been screened and annotated using
\code{LOBSTAHS} function \code{\link{doLOBscreen}}.
}
\section{Objects from the class}{
Objects can be created with the \code{\link{LOBSet}} constructor. A 'LOBSet' can
also be created by calls of the form \code{new("LOBSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{peakdata}:}{Object of class \code{"data.frame"}, containing
peakdata by compound assignment }
    \item{\code{iso_C3r}:}{Object of class \code{"list"}, a list of the
\code{match_ID}s of possible regioisomers of each compound }
    \item{\code{iso_C3f}:}{Object of class \code{"list"}, a list of the
\code{match_ID}s of possible functional structural isomers of each compound }
    \item{\code{iso_C3c}:}{Object of class \code{"list"}, a list of the
\code{match_ID}s of possible isobars of this compound }
    \item{\code{LOBscreen_diagnostics}:}{Object of class \code{"data.frame"},
numbers of peaks, peakgroups, adducts, and unique parent compounds present in
the dataset after application of each screening criterion in
\code{\link{doLOBscreen}} }
    \item{\code{LOBisoID_diagnostics}:}{Object of class \code{"data.frame"},
numbers of peakgroups and parent compounds to which the various isomer
annotations have been assigned }
    \item{\code{LOBscreen_settings}:}{Object of class \code{"list"}, captures
the settings used in \code{\link{doLOBscreen}} to generate the \code{"LOBSet"} }
    \item{\code{polarity}:}{Object of class \code{"factor"}, polarity of data in
the \code{"LOBSet"} }
    \item{\code{sampnames}:}{Object of class \code{"character"}, the names of
the samples from which the \code{"LOBSet"} was generated }
  }
}
\section{Methods}{
  \describe{
    \item{show}{
      \code{signature(object = "LOBSet")}: ... 
    }
    \item{LOBisoID_diagnostics}{
      \code{signature(object = "LOBSet")}: get \code{LOBisoID_diagnostics} slot
    }
    \item{LOBscreen_diagnostics}{
      \code{signature(object = "LOBSet")}: get \code{LOBscreen_diagnostics} slot
    }
    \item{LOBscreen_settings}{
      \code{signature(object = "LOBSet")}: get \code{LOBscreen_settings} slot
    }
    \item{sampnames}{
      \code{signature(object = "LOBSet")}: get \code{sampnames} slot
    }
    \item{peakdata}{
      \code{signature(object = "LOBSet")}: get \code{peakdata} slot
    }
    \item{polarity}{
      \code{signature(object = "LOBSet")}: get \code{polarity} slot
    }
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
  \code{\link{doLOBscreen}},
  \code{\link{getLOBpeaklist}},
  \code{\link{LOBSet}},
  \code{\linkS4class{xsAnnotate}}
}
\examples{
## return object information

showClass("LOBSet")
}
\keyword{classes}