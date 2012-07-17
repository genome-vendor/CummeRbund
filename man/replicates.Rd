\name{replicates}
\alias{replicates}
\alias{replicates,CuffSet-method}
\alias{replicates,CuffData-method}
\alias{replicates,CuffFeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get replicate sample list from CuffData object
}
\description{
Returns a list of replicate names from a CuffData or CuffFeatureSet object
}
\usage{
\S4method{replicates}{CuffData}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffSet','CuffData')
}
}
\details{
None
}
\value{
A list of replicate sample names
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
	replicates(a@genes)
}

