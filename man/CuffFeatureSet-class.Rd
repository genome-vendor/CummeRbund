\name{CuffFeatureSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffFeatureSet-class}
\alias{diffData,CuffFeatureSet-method}
\alias{annotation}
\alias{annotation,CuffFeatureSet-method}
\alias{featureNames,CuffFeatureSet-method}
\alias{features,CuffFeatureSet-method}
\alias{fpkmMatrix,CuffFeatureSet-method}
\alias{repFpkmMatrix,CuffFeatureSet-method}
\alias{countMatrix,CuffFeatureSet-method}
\alias{samples,CuffFeatureSet-method}
\alias{length,CuffFeatureSet-method}

\title{Class "CuffFeatureSet"}
\description{
A 'data' container class for all FPKM, annotation, and differential expression data for a set of features (genes, isoforms, TSS, CDS).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffFeatureSet", annotation, fpkm, diff, ...)}.

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
    \item{diffData}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{featureNames}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{features}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{fpkmMatrix}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{repFpkmMatrix}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{countMatrix}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{samples}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{length}{\code{signature(object = "CuffFeatureSet")}: ... }
	 }
}
\section{Accessors}{
	\describe{
		\item{annotation}{\code{signature(object="CuffFeatureSet")}: Access @annotation slot}
	}
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}



\seealso{
\code{\linkS4class{CuffGeneSet}}

}
\examples{
showClass("CuffFeatureSet")
}
\keyword{classes}
