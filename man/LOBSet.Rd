\name{LOBSet}
\Rdversion{1.2}
\alias{LOBSet}

\title{LOBSet constructor for manual creation or reconstruction of a LOBSet
object}
\description{
Constructor function for manual creation or reconstruction of a
\code{LOBSet} object for HPLC-MS peak data that have been screened using
\code{LOBSTAHS}.
}
\usage{
LOBSet(peakdata = NULL, iso_C3r = NULL, iso_C3f = NULL, iso_C3c = NULL,
       LOBscreen_diagnostics = NULL, LOBisoID_diagnostics = NULL,
       LOBscreen_settings = NULL, polarity = c("positive","negative"),
       sampnames = NULL)
}
\arguments{
  \item{peakdata}{
An object of class \code{"data.frame"} containing peak data and LOBSTAHS
annotation information by assignment. Column headings and data types should
conform to those of the peakdata slot of a \code{\link{LOBSet-class}} object
produced using the function \code{\link{doLOBscreen}}. Each row in the
peakdata table represents one compound assignment made by LOBSTAHS. The format
of the data frame occupying the peakdata slot in a \code{LOBSet} object can be
obtained using the \code{\link{peakdata}} accessor for objects of
\code{\link{LOBSet-class}}.
}
  \item{iso_C3r}{
An object of class \code{"list"}, a list of the \code{match_ID}s of possible
regioisomers of each compound for which there is a row in \code{peakdata}.
Length of \code{iso_C3r} should equal the number of rows in the data frame
given for \code{peakdata}.
}
  \item{iso_C3f}{
An object of class \code{"list"}, a list of the \code{match_ID}s of possible
functional structural isomers of each compound for which there is a row in
\code{peakdata}. Length of \code{iso_C3f} should equal the number of rows in the
data frame given for \code{peakdata}.
}
  \item{iso_C3c}{
An object of class \code{"list"}, a list of the \code{match_ID}s of the isobars
of each compound for which there is a row in \code{peakdata}. Length of
\code{iso_C3c} should equal the number of rows in the data frame given for
\code{peakdata}.
}
  \item{LOBscreen_diagnostics}{
An object of class \code{"data.frame"}, containing diagnostic information
recorded by the function \code{\link{doLOBscreen}} during creation of a
\code{\link{LOBSet}}. The numbers of peaks, peakgroups, adducts, and unique
parent compounds present the dataset after application of each LOBSTAHS
screening criterion.
}
  \item{LOBisoID_diagnostics}{
An object of class \code{"data.frame"}, containing isomer and isobar summary
statistics. The numbers of peakgroups and parent compounds to which the various
isomer identifications have been applied by \code{\link{doLOBscreen}} to a given
\code{LOBSet}.
}
  \item{LOBscreen_settings}{
An object of class \code{"list"}; the settings used in \code{\link{doLOBscreen}}
to generate the \code{LOBSet}
}
  \item{polarity}{
An object of class \code{"factor"}, polarity of data in the \code{LOBSet}. Must
be either "negative" or "positive."
}
  \item{sampnames}{
An object of class \code{"character"}; the names of the samples from which the
\code{LOBSet} was generated.
}
}
\details{Typically, a \code{LOBSet} will be created from a CAMERA
\code{\link{xsAnnotate-class}} object using the \code{LOBSTAHS} function
\code{\link{doLOBscreen}}. The \code{LOBSet} constructor function is therefore
provided only for manual object creation; it will not be needed by most users.
}
\value{
  A \code{"LOBSet-class"} object.
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
  \code{\link{doLOBscreen}},
  \code{\linkS4class{xsAnnotate}}
}
\examples{
## create an empty LOBSet for positive ion mode data

myLOBSet = LOBSet(peakdata = NULL, iso_C3r = NULL, iso_C3f = NULL,
                  iso_C3c = NULL, LOBscreen_diagnostics = NULL,
                  LOBisoID_diagnostics = NULL, LOBscreen_settings = NULL,
                  polarity = "positive", sampnames = NULL)
}
\keyword{methods}
