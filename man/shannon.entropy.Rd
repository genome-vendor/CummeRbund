\name{shannon.entropy}
\alias{shannon.entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shannon entropy
}
\description{
Calculates the Shannon entropy for a probability distribution
}
\usage{
shannon.entropy(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of probabilities (must sum to ~1)
}
}
\details{
None
}
\value{
Returns a numeric value for the Shannon entropy of the supplied probability distribution
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

\examples{
	x<-sample(1:500,50)
	p<-x/sum(x)
	shannon.entropy(p)
}
