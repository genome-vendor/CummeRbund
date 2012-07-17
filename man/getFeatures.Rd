\name{getFeatures}
\alias{getFeatures}
\alias{getFeatures,CuffSet-method}

\title{
getFeatures
}
\description{
Primary accessor from a CuffSet object to retrieve all related information for >1 (MANY) given FEATURES, indexed by tracking id.
}
\usage{
\S4method{getFeatures}{CuffSet}(object, featureIdList, sampleIdList=NULL,level='isoforms')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class 'CuffSet' (Primary 'pointer' object for Cufflinks data).
}
  \item{featureIdList}{
A vector of 'isoform_id', 'TSS_group_id', or 'CDS_id' to identify which features for which you would like to retrieve all information.
}
  \item{sampleIdList}{
A vector of sample names used to subset or re-order samples in returned object
}
  \item{level}{
Feature level to be queried for significance (must be one of c('isoforms','TSS','CDS')
}
}
\details{
None.
}
\value{
Returns a CuffFeatureSet object containing all related information for a given set of tracking_id values
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
	sample.isoform.ids<-sample(featureNames(isoforms(cuff)),10)
	myGene<-getFeatures(cuff,sample.isoform.ids) # Retrieve all information for a set of 10 sampled features.
}