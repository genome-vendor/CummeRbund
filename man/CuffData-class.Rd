\name{CuffData-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffData-class}
\alias{dim,CuffData-method}
\alias{getFeatures,CuffData-method}

\title{Class "CuffData"}
\description{
A 'pointer' class for all information (FPKM, annotation, differential expression) for a given feature type (genes, isoforms, TSS, CDS). The methods for this function
communicate directly with the SQL backend to present data to the user.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffData", DB, tables, filters, type, idField, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{tables}:}{Object of class \code{"list"} ~~ }
    \item{\code{filters}:}{Object of class \code{"list"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{idField}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "CuffData")}: ... }
    \item{getFeatures}{\code{signature(object = "CuffData")}: ... }
	 }
}
\references{
None
}
\author{
Loyal A. Goff
}
\note{
None
}



\seealso{
None
}
\examples{
showClass("CuffData")
}
\keyword{classes}
