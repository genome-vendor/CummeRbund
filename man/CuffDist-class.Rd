\name{CuffDist-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffDist-class}
\alias{dim,CuffDist-method}
\alias{samples,CuffDist-method}

\title{Class "CuffDist"}
\description{
A 'pointer' class to information relative to the distribution-level tests (promoters, splicing, and relative CDS usage)
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffDist", DB, table, type, idField, ...)}.
 
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{table}:}{Object of class \code{"character"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{idField}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "CuffDist")}: ... }
    \item{samples}{\code{signature(x = "CuffDist")}: ... }
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
showClass("CuffDist")
}
\keyword{classes}
