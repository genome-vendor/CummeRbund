\name{csScatter}
\alias{csScatter}
\alias{csScatter,CuffData-method}
\alias{csScatter,CuffFeatureSet-method}

\title{
Scatter Plot
}
\description{
A scatter plot comparing the FPKM  values from two samples in a cuffdiff run.
}
\usage{
\S4method{csScatter}{CuffData}(object, x, y, logMode=TRUE, pseudocount=1.0, labels, smooth=FALSE, colorByStatus = FALSE, drawRug=TRUE, ...)
}

\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet')
}
  \item{x}{
Sample name for x axis
}
  \item{y}{
Sample name for y axis
}
  \item{logMode}{
Logical argument to log2-transform data (default: T )
}
  \item{pseudocount}{
Value to add to zero FPKM values for log transformation (default: 0.0001)
}
  \item{smooth}{
Logical argument to add a smooth-fit regression line
}
  \item{labels}{
A list of tracking_ids or gene_short_names that will be 'callout' points in the plot for reference. Useful for finding genes of interest in the field. Not implemented yet.
}
  \item{colorByStatus}{
A logical argument whether or not to color the points by 'significant' Y or N. [Default = FALSE]
}
  \item{drawRug}{
A logical argument whether or not to draw the rug for x and y axes [Default = TRUE]
}
  \item{\dots}{
Additional arguments to csScatter
}
}
\details{
None
}
\value{
ggplot object with geom_point and geom_rug layers
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
	genes<-a@genes #Create CuffData object for all genes
	s<-csScatter(genes,'hESC','Fibroblasts',smooth=TRUE) #Create plot object
	s #render plot object
}
