\name{csSpecificity}
\alias{csSpecificity}
\alias{csSpecificity,CuffFeatureSet-method}
\alias{csSpecificity,CuffData-method}

\title{
csSpecificity
}
\description{
Returns a matrix of 'Specificity scores' (S) defined as 1-JSD(p_g,q_i) where p_g is the Log10+1 expression profile of a gene g across all conditions j, collapsed into a probability distribution,
and q_i is the unit vector of 'perfect expression' in a given condition i.
}
\usage{
\S4method{csSpecificity}{CuffFeatureSet}(object,logMode=T,pseudocount=1,relative=FALSE,...)
\S4method{csSpecificity}{CuffData}(object,logMode=T,pseudocount=1,relative=FALSE,...)
}
\arguments{
  \item{object}{
An object of class CuffFeatureSet, CuffGeneSet, or CuffData.
}
  \item{logMode}{
A logical argument to log10-transform FPKM values prior to plotting.
}
  \item{pseudocount}{
Value to be added to FPKM for appropriate log transformation and clustering. (Avoids zero-based errors)
}
  \item{relative}{
A logical argument that when TRUE, will scale the S values from 0-1 by dividing by max(S)
}
  \item{\dots}{
Additional arguments to fpkmMatrix.
}
}
\details{
None
}
\author{
Loyal A. Goff
}
\source{
None
}
\references{
None.
}
\examples{
	data(sampleData)
	csSpecificity(sampleGeneSet)
}
\keyword{datasets}
