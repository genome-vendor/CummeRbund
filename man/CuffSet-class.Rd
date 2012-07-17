\name{CuffSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffSet-class}
%Accessors
\alias{genes}
\alias{isoforms}
\alias{TSS}
\alias{CDS}
\alias{promoters}
\alias{splicing}
\alias{relCDS}

\alias{genes,CuffSet-method}
\alias{isoforms,CuffSet-method}
\alias{TSS,CuffSet-method}
\alias{CDS,CuffSet-method}
\alias{promoters,CuffSet-method}
\alias{splicing,CuffSet-method}
\alias{relCDS,CuffSet-method}

\alias{[,CuffSet-method}

\title{Class "CuffSet"}
\description{
A 'pointer' class to connect to, and retrieve data from the SQLite backend database. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffSet", DB, conditions, genes, isoforms, TSS, CDS, promoters, splicing, relCDS, ...)}.
Available methods are primary accessors to retrieve CuffGeneSet or CuffGene objects for manipulation.
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{conditions}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{genes}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{promoters}:}{Object of class \code{"CuffDist"} ~~ }
    \item{\code{runInfo}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{splicing}:}{Object of class \code{"CuffDist"} ~~ }
    \item{\code{relCDS}:}{Object of class \code{"CuffDist"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "CuffSet")}: ... }
	 }
}

\section{Accessors}{
	\describe{
		\item{genes}{\code{signature(object = "CuffSet")}: Access @genes slot }
		\item{isoforms}{\code{signature(object = "CuffSet")}: Access @isoforms slot }
		\item{TSS}{\code{signature(object = "CuffSet")}: Access @TSS slot }
		\item{CDS}{\code{signature(object = "CuffSet")}: Access @CDS slot }
		\item{promoters}{\code{signature(object = "CuffSet")}: Access @promoters slot }
		\item{splicing}{\code{signature(object = "CuffSet")}: Access @splicing slot }
		\item{relCDS}{\code{signature(object = "CuffSet")}: Access @relCDS slot }
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
None.
}
\examples{
showClass("CuffSet")
}
\keyword{classes}
