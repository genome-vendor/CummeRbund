\name{CuffFeature-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffFeature-class}
\alias{annotation,CuffFeature-method}
\alias{diffData,CuffFeature-method}
\alias{fpkmMatrix,CuffFeature-method}
\alias{length,CuffFeature-method}
\alias{samples,CuffFeature-method}

\title{Class "CuffFeature"}
\description{
A 'data' container class for all FPKM, annotation, and differential expression data for a single feature (gene, isoform, TSS, or CDS).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffFeature", annotation, fpkm, diff, ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{repFpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{count}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{fpkmMatrix}{\code{signature(object="CuffFeature")}: ...}
    \item{length}{\code{signature(x = "CuffFeature")}: ... }
	 }
}

\section{Accessors}{
	\describe{
		\item{annotation}{\code{signature(object="CuffFeature")}: Access @annotation slot}
		\item{diffData}{\code{signature(object="CuffFeature")}: Access @diff slot}
		\item{samples}{\code{signature(object="CuffFeature")}: Get vector of samples}
	}
}
\references{
None
}
\author{
Loyal A. Goff
}
\note{
'CuffGene' is a superclass of 'CuffFeature' that links gene information for a given gene with all isoform-, TSS-, and CDS-level data for the given gene.
}



\seealso{
\code{\linkS4class{CuffGene}}
}
\examples{
showClass("CuffFeature")
}
\keyword{classes}
