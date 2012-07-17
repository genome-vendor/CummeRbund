\name{featureNames}
\alias{featureNames}
\alias{featureNames,CuffData-method}

\title{
Feature names
}
\description{
Retrive a vector of feature names from a 'CuffData' or 'CuffFeatureSet' object
}
\usage{
\S4method{featureNames}{CuffData}(object)
}

\arguments{
  \item{object}{
An object of class ('CuffData' or 'CuffFeatureSet')
}
}
\details{
None
}
\value{
A list of feature names
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
featureNames(sampleGeneSet)
}
