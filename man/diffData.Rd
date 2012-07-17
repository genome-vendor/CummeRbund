\name{diffData}
\alias{diffData}
\alias{diffData,CuffData-method}

\title{
Differential comparison data
}
\description{
An accessor method to retrieve differential expression data from a 'CuffData', 'CuffFeatureSet', or 'CuffFeature' object
}
\usage{
\S4method{diffData}{CuffData}(object, x, y, features=FALSE)
}

\arguments{
  \item{object}{
An object of class ('CuffData' or 'CuffFeatureSet')
}
  \item{x}{
  Optional, if x and y are both missing, data for all pairwise differential comparisons are returned, otherwise if x and y are sample names from the 'samples' table,
  than only differential data pertaining to those two samples are returned.
}
  \item{y}{
See 'x'
}
  \item{features}{
A logical value that returns all feature-level data as part of data.frame when true.  object must be of class 'CuffData'.
}

  \item{\dots}{
Additional arguments.
}
}
\details{
None
}
\value{
A data.frame object
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
	data(sampleData)
	diff<-diffData(sampleGeneSet) #returns a dataframe of differential expression data from sample CuffGeneSet object.
}

