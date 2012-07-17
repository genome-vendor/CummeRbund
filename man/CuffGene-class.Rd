\name{CuffGene-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffGene-class}
\alias{genes,CuffGene-method}
\alias{isoforms,CuffGene-method}
\alias{TSS,CuffGene-method}
\alias{CDS,CuffGene-method}
\alias{promoters,CuffGene-method}
\alias{relCDS,CuffGene-method}
\alias{splicing,CuffGene-method}
\alias{length,CuffGene-method}

\title{Class "CuffGene"}
\description{
A 'data' container class for all FPKM, annotation, and differential expression Data (as well as for all linked features) for a given gene.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffGene", id, isoforms, TSS, CDS, promoters, splicing, relCDS, annotation, fpkm, diff, ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{promoters}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{relCDS}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{splicing}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CuffFeature}"}, directly.
}
\section{Methods}{
\describe{
	\item{length}{\code{signature(object="CuffFeature")}: Part of length validation (internal use only)}
}
}

\section{Accessors}{
	\describe{
		\item{genes}{\code{signature(object = "CuffGene")}: Access @genes slot }
		\item{isoforms}{\code{signature(object = "CuffGene")}: Access @isoforms slot }
		\item{TSS}{\code{signature(object = "CuffGene")}: Access @TSS slot }
		\item{CDS}{\code{signature(object = "CuffGene")}: Access @CDS slot }
		\item{promoters}{\code{signature(object = "CuffGene")}: Access @CDS slot }
		\item{relCDS}{\code{signature(object = "CuffGene")}: Access @CDS slot }
		\item{splicing}{\code{signature(object = "CuffGene")}: Access @CDS slot }
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
\code{\linkS4class{CuffFeature}}
}
\examples{
showClass("CuffGene")
}
\keyword{classes}
