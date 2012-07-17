\name{distValues}
\alias{distValues}
\alias{distValues,CuffDist-method}

\title{
distValues
}
\description{
Returns a data.frame of distribution-level test values from a CuffDist object (@promoters, @splicing, @relCDS)
}
\usage{
\S4method{distValues}{CuffDist}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class 'CuffDist'
}
  \item{\dots}{
Additional arguments to distValues
}
}
\details{
None
}
\value{
Returns a data.frame of distribution-level test values.
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) # Read cufflinks data and create CuffSet object
	distValues(a@promoters) # returns data.frame of values from CuffDist object in slot 'promoters'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

