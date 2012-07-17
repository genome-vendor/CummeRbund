\name{csBoxplot}
\alias{csBoxplot}
\alias{csBoxplot,CuffData-method}

\title{
csBoxplot
}
\description{
Creates a ggplot2 plot object with a geom_box layer displaying summary statistics for FPKM values across samples (x).
}
\usage{
\S4method{csBoxplot}{CuffData}(object, logMode=TRUE, pseudocount=0.0001, replicates=FALSE,...)
}

\arguments{
  \item{object}{
An object of class CuffData.
}
  \item{logMode}{
A logical argument to log10 -transform FPKM values.
}
  \item{pseudocount}{
Value added to FPKM to avoid log-transform issues. 
}
  \item{replicates}{
A logical value whether or not to plot individual replicates or aggregate condition values.
}

  \item{\dots}{
Additional arguments to csBoxplot
}
}
\details{
None
}
\value{
A ggplot2 plot object with a geom_box layer.
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create CuffSet object
	genes<-a@genes #CuffData object for all genes
	csBoxplot(genes)
}

