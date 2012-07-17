\name{runInfo}
\alias{runInfo}
\alias{runInfo,CuffSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve run parameters and information from a CuffSet object
}
\description{
Returns a data.frame of cuffdiff run parameters and information
}
\usage{
\S4method{runInfo}{CuffSet}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffSet')
}
}
\details{
None
}
\value{
A data.frame of run parameters
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
None
}
\examples{

	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object
	runInfo(a)
}

