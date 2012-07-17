\name{getGeneId}
\alias{getGeneId}
\alias{getGeneId,CuffSet-method}

\title{
getGene
}
\description{
A helper function to retrieve the gene_ids for a given list of feature ids (e.g. isoform_ids, tss_group_ids, or CDS_ids). This should not be
called directly by the user
}
\usage{
\S4method{getGeneId}{CuffSet}(object, idList)
}

\arguments{
  \item{object}{
An object of class 'CuffSet' (Primary 'pointer' object for Cufflinks data).
}
  \item{idList}{
A character string to identify the identifiers for which you would like to retrieve corresponding gene_id values.
}
}
\details{
None.
}
\value{
Returns a vector of gene_id values corresponding to genes from which idList are sub-features.
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
	cuff<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create master CuffSet object
	sampleFeatureIds<-sample(featureNames(isoforms(cuff)),10)
	correspondingGeneIds<-getGeneId(cuff,sampleFeatureIds) # Retrieve gene_id values for parent genes of sampleFeatureIds.
}
