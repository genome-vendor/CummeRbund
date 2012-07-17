\name{MAplot}
\alias{MAplot}
\alias{MAplot,CuffData-method}

\title{
MAplot
}
\description{
Creates an M vs A plot (Avg intensity vs log ratio) for a given pair of conditions across all fpkms
}
\usage{
\S4method{MAplot}{CuffData}(object,x,y,logMode=T,pseudocount=1,smooth=FALSE,useCount=FALSE)
}

\arguments{
  \item{object}{
An object of class 'CuffData'.
}
  \item{x}{
Sample name from 'samples' table for comparison
}
  \item{y}{
Sample name from 'samples' table for comparison
}
  \item{logMode}{
A logical argument to log10-transform FPKM values prior to plotting.
}
  \item{pseudocount}{
Value to be added to FPKM for appropriate log transformation and clustering. (Avoids zero-based errors)
}
  \item{smooth}{
Logical argument whether or not to draw a smoothed line fit through data.
}
  \item{useCount}{
Logical argument whether or not to use mean counts instead of FPKM values.
}
 
 
}
\details{
None
}
\value{
Returns a ggplot MvsA plot object.
}
\references{
None.
}
\author{
Loyal A. Goff and Cole Trapnell
}
\note{
None
}

\examples{
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object from sample data
	genes<-a@genes #Create CuffData object for all 'genes'
	d<-MAplot(genes,'hESC','Fibroblasts') #Create csDensity plot
	d #Render plot
}

\keyword{heatmap}
