\name{findSimilar}
\alias{findSimilar}
\alias{findSimilar,CuffSet-method}

\title{
findSimilar
}
\description{
Returns a CuffGeneSet containing n genes with the most similar expression profiles to gene/profile x.
}
\usage{
\S4method{findSimilar}{CuffSet}(object, x, n,distThresh,returnGeneSet=TRUE,...)
}

\arguments{
  \item{object}{
A object of class 'CuffSet'
}
  \item{x}{
A 'gene_id' or 'gene_short_name' from which to look up an expression profile OR a vector of expression values to compare all genes (vector must have same length and order of 'samples')
}
  \item{n}{
Number of similar genes to return
}
  \item{distThresh}{
A thresholding value on which to filter results based on JS-distance (e.g. A distThresh of 1.0 will return all genes, 0.0 will return those genes with 'perfect identity' to the gene of interest.) 
}
  \item{returnGeneSet}{
A logical value whether to return a CuffGeneSet object [default] or a distance-ranked data frame of similar genes.  The latter is useful if you want to explore the returned list based on distances.
}
  \item{\dots}{
Additional arguments to fpkmMatrix call within findSimilar (e.g. fullnames=T)
}
}
\details{
By default, returns a CuffGeneSet object with n similar genes.  This may change in the future.
}
\value{
A CuffGeneSet object of n most similar genes to x.
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund"))
	mySimilarGenes<-findSimilar(a,"PINK1",10)
}