\name{csClusterPlot}
\alias{csClusterPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
csClusterPlot
}
\description{
Replaces the default plotting behavior of the old csCluster.  Takes as an argument the output of csCluster and plots expression profiles of features facet by cluster.
}
\usage{
csClusterPlot(clustering, pseudocount=1.0,drawSummary=TRUE,sumFun=mean_cl_boot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
The output of csCluster. (Must be the output of csCluster. Only this data format contains the necessary information for csClusterPlot.)
}
  \item{pseudocount}{
Value added to FPKM to avoid log transformation issues.
}
  \item{drawSummary}{
Logical value whether or not to draw a summary line for each cluster (by default this is the cluster mean)
}
  \item{sumFun}{
Summary function used to by drawSummary (default: mean_cl_boot)
}
}
\details{
This replaces the default plotting behavior of the old csCluster() method. This was necessary so as to preserve the cluster information obtained by csCluster in a stable format.
The output of csClusterPlot is a ggplot2 object of expressionProfiles faceted by cluster ID.

}		
\value{
A ggplot2 object of expressionProfiles faceted by cluster ID.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}


\examples{
	data(sampleData)
	myClustering<-csCluster(sampleGeneSet,k=4)
	csClusterPlot(myClustering)
}

