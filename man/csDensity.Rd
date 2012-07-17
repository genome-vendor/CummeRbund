\name{csDensity}
\alias{csDensity}
\alias{csDensity,CuffData-method}
\alias{csDensity,CuffFeatureSet-method}

\title{
Density plot of CuffData
}
\description{
Creates a smoothed density plot, by sample, for log10 FPKM values from a cuffdiff run.
}
\usage{
\S4method{csDensity}{CuffData}(object, logMode=TRUE, pseudocount=1.0, labels, features=FALSE, replicates=FALSE,...)
\S4method{csDensity}{CuffFeatureSet}(object, logMode=TRUE, pseudocount=1.0, labels, features=FALSE, replicates=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class CuffData.
}
  \item{logMode}{
A logical value of whether or not to log10-transform FPKM values.  By default this is TRUE.
}
  \item{pseudocount}{
Pseudocount value added to FPKM to avoid errors in log-transformation of true zero values.
}
  \item{labels}{
A list of tracking_id values or gene_short_name values used for 'callout' points on the density plot for reference. (Not implemented yet).
}
  \item{features}{
Will include all fields from 'features' slot in returned ggplot object.
Useful for further manipulations of plot object using feature-level attributes (e.g. gene_type, class_code, etc)
}
  \item{replicates}{
A logical value whether or not to plot individual replicates or aggregate condition values.
}
  \item{\dots}{
Additional arguments
}
}
\details{
Creates a density plot, by sample, for log10-transformed FPKM values from a cuffdiff run.
}
\value{
A ggplot2 plot object
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object from sample data
	genes<-a@genes #Create CuffData object for all 'genes'
	d<-csDensity(genes) #Create csDensity plot
	d #Render plot
}
