\name{makeprobsvec}
\alias{makeprobsvec}

\title{
makeprobsvec
}
\description{
Sums a vector of numerics and divides by the sum
}
\usage{
makeprobsvec(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of numerics
}
}
\details{
None
}
\value{
A vector of probabilities
}
\references{
None
}
\author{
Loyal A. Goff
}
\note{
Should not be called directly by user.
}
\examples{
	p<-sample(1:5000,20)
	makeprobsvec(p)
}
