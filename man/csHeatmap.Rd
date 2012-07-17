\name{csHeatmap}
\alias{csHeatmap}
\alias{csHeatmap,CuffFeatureSet-method}

\title{
csHeatmap
}
\description{
Creates a ggplot plot object with a geom_tile layer of FPKM values per feature and sample.
}
\usage{
\S4method{csHeatmap}{CuffFeatureSet}(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, 
		border=FALSE, heatscale= c(low='lightyellow',mid='orange',high='darkred'), heatMidpoint=NULL, fullnames = T, replicates=FALSE, ...)
}

\arguments{
  \item{object}{
An object of class 'CuffFeatureSet' or 'CuffGeneSet'
}
  \item{rescaling}{
Rescaling can either be 'row' or 'column' OR you can pass rescale a function that operates on a matrix to do your own rescaling. Default is 'none'.
}
  \item{clustering}{
Clustering can either be 'row','column','none', or 'both', in which case the appropriate indices are re-ordered based on the pairwise Jensen-Shannon distance of FPKM values.
Alternatively you can pass your own clustering function so long as the returned value is a re-ordered matrix.
}
  \item{labCol}{
A logical argument to display column labels.
}
  \item{labRow}{
A logical argument to display row labels.
}
  \item{logMode}{
A logical argument to log10-transform FPKM values prior to plotting.
}
  \item{pseudocount}{
Value to be added to FPKM for appropriate log transformation and clustering. (Avoids zero-based errors)
}
  \item{border}{
A logical argument to draw border around plot.
}
  \item{heatscale}{
A list with min length=2, max length=3 that detail the low,mid,and high colors to build the color scale.
}
  \item{heatMidpoint}{
Value for midpoint of color scale.
}
  \item{fullnames}{
A logical value whether to use 'fullnames' (concatenated gene_short_name and gene_id) for rows in heatmap. Default [ TRUE ].
}
  \item{replicates}{
A logical value whether or not to plot individual replicates or aggregate condition values.
}
  \item{\dots}{
Additional arguments to csHeatmap
}
}
\details{
None
}
\value{
A ggplot2 plot object with a geom_tile layer to display FPKM values by sample (x) and feature (y)
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
	data(sampleData)
	csHeatmap(sampleGeneSet)
}

\keyword{heatmap}
