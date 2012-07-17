\name{JSdistVec}
\alias{JSdistVec}

\title{
JSdistVec
}
\description{
Returns the Jensen-Shannon Distance (square root of JS divergence) between two probability vectors.
}
\usage{
JSdistVec(p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of probabilities
}
  \item{q}{
A vector of probabilities
}
}
\details{
Should not be called directly by user.
}
\value{
Returns the JS distance as a numeric
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
	p<-sample(1:5000,20)
	q<-sample(1:5000,20)
	p<-makeprobsvec(p)
	q<-makeprobsvec(q)
	JSdistVec(p,q)
	
}