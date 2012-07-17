\name{repFpkmMatrix}
\alias{repFpkmMatrix}
\alias{repFpkmMatrix,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
repFpkmMatrix
}
\description{
Retrieve repFpkm values as gene by condition matrix
}
\usage{
\S4method{repFpkmMatrix}{CuffData}(object,fullnames=FALSE,repIdList)
}

\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet','CuffGeneSet','CuffGene',or 'CuffFeature')
}
  \item{fullnames}{
A logical value whether or not to concatenate gene_short_name and tracking_id values (easier to read labels)
}
  \item{repIdList}{
A vector of sample names to subset the resulting matrix.
}
}
\details{
None.
}
\value{
A feature x condition matrix of replicate FPKM values.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None
}

\examples{
	data(sampleData)
	repFpkmMatrix(sampleGeneSet)
}
