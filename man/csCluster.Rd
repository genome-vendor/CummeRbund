\name{csCluster}
\alias{csCluster}
\alias{csCluster,CuffFeatureSet-method}

\title{
csCluster
}
\description{
Returns a ggplot2 plot object with geom_line layer plotting FPKM values over conditions faceted by k-means clustering clusters. (Euclidean).  This is very crude at this point.
This does not return any of the clustering information directly, but if you want it, you can retrieve it from the ggplot object returned.
}
\usage{
\S4method{csCluster}{CuffFeatureSet}(object,k,pseudocount=1,...)
}
\arguments{
  \item{object}{
An object of class CuffFeatureSet.
}
  \item{k}{
Number of pre-defined clusters to attempt to find.
}
  \item{pseudocount}{
Value added to FPKM to avoid log-transform issues.
}
  \item{\dots}{
Additional arguments to pam.
}
}
\details{
Uses 'kmeans' function.
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
	csCluster(sampleGeneSet,4)
}
\keyword{datasets}
