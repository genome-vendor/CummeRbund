\name{dispersionPlot}
\alias{dispersionPlot}
\alias{dispersionPlot,CuffData-method}

\title{
Mean count vs dispersion plot
}
\description{
A scatter plot comparing the mean counts against the estimated dispersion for a given level of features from a cuffdiff run.
}
\usage{
\S4method{dispersionPlot}{CuffData}(object)
}

\arguments{
  \item{object}{
An object of class ('CuffData')
}

}
\details{
None
}
\value{
ggplot object with geom_point layer
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
	genes<-genes(a) #Create CuffData object for all genes
	d<-dispersionPlot(genes) #Create plot object
	d #render plot object
}
