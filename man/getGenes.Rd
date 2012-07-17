\name{getGenes}
\alias{getGenes}
\alias{getGenes,CuffSet-method}

\title{
getGenes
}
\description{
Primary accessor from a CuffSet object to retrive all related information for >1 (MANY) given genes, indexed by gene_id or gene_short_name.
}
\usage{
\S4method{getGenes}{CuffSet}(object, geneIdList, sampleIdList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class 'CuffSet' (Primary 'pointer' object for Cufflinks data).
}
  \item{geneIdList}{
A vector of gene_ids or gene_short_namesto identify which genes for which you would like to retrieve all information.
}
  \item{sampleIdList}{
A vector of sample names used to subset or re-order samples in returned object
}
}
\details{
None.
}
\value{
Returns a CuffGeneSet object containing all related information for a given set of gene_id or gene_short_name values
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
Right now, this does not return an error if it cannot find a gene.  (this is probably a bad thing...)
}

\examples{
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create master CuffSet object
	data(sampleData)
	sampleIDs
	myGene<-getGenes(a,sampleIDs) # Retrieve all information for a set of 20 'sample' genes.
}