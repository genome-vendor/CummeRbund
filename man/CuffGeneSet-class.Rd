\name{CuffGeneSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffGeneSet-class}
\alias{genes,CuffGeneSet-method}
\alias{isoforms,CuffGeneSet-method}
\alias{TSS,CuffGeneSet-method}
\alias{CDS,CuffGeneSet-method}
\alias{promoters,CuffGeneSet-method}
\alias{relCDS,CuffGeneSet-method}
\alias{splicing,CuffGeneSet-method}

\title{Class "CuffGeneSet"}
\description{
A 'data' container class for all FPKM, annotation, and differential expression data (an associated features) for a given set of genes.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffGeneSet", annotation, fpkm, diff, ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{ids}:}{Object of class \code{"character"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{promoters}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{relCDS}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{splicing}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CuffFeatureSet}"}, directly.
}
\section{Methods}{
No methods defined with class "CuffGeneSet" in the signature.
}

\section{Accessors}{
	\describe{
		\item{genes}{\code{signature(object = "CuffGeneSet")}: Access @genes slot }
		\item{isoforms}{\code{signature(object = "CuffGeneSet")}: Access @isoforms slot }
		\item{TSS}{\code{signature(object = "CuffGeneSet")}: Access @TSS slot }
		\item{CDS}{\code{signature(object = "CuffGeneSet")}: Access @CDS slot }
		\item{promoters}{\code{signature(object = "CuffGeneSet")}: Access @promoters slot }
		\item{relCDS}{\code{signature(object = "CuffGeneSet")}: Access @relCDS slot }
		\item{splicing}{\code{signature(object = "CuffGeneSet")}: Access @splicing slot }
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
\code{\linkS4class{CuffFeatureSet}}
}
\examples{
showClass("CuffGeneSet")
}
\keyword{classes}
