\name{features}
\alias{features}
\alias{features,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Features
}
\description{
Returns a data frame of features from a CuffData, CuffFeatureSet, or CuffFeature object
}
\usage{
\S4method{features}{CuffData}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet', or 'CuffFeature')
}
}
\details{
If features have been added to cuffData via addFeatures() then these will be presented as well.
}
\value{
A data.frame of feature-level information
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
	features(sampleGeneSet)
}
