\name{expressionBarplot}
\alias{expressionBarplot}
\alias{expressionBarplot,CuffFeature-method}
\alias{expressionBarplot,CuffFeatureSet-method}

\title{
Barplot
}
\description{
A barplot of FPKM values with confidence intervals for a given gene, set of genes, or features of a gene (e.g. isoforms, TSS, CDS, etc).
}
\usage{
\S4method{expressionBarplot}{CuffFeatureSet}(object, logMode=TRUE, pseudocount=1.0, showErrorbars=TRUE, showStatus=TRUE, replicates=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffFeatureSet','CuffGeneSet','CuffFeature','CuffGene')
}
  \item{logMode}{
A logical value whether or not to draw y-axis on log10 scale. Default = FALSE.
}
  \item{pseudocount}{
Numerical value added to each FPKM during log-transformation to avoid errors.
}
  \item{showErrorbars}{
A logical value whether or not to draw error bars. Default = TRUE
}
  \item{showStatus}{
A logical value whether or not to draw visual queues for quantification status of a given gene:condition. Default = TRUE
}
  \item{replicates}{
A logical value whether or not to plot individual replicates or aggregate condition values.
}
  \item{\dots}{
Additional arguments.
}
}
\details{
None
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
Need to implement logMode and features for this plotting method.
}

\examples{
	data(sampleData)
	PINK1 # sample CuffFeature object
	expressionBarplot(PINK1) #Barplot of PINK1 FPKM values
	expressionBarplot(PINK1@isoforms) #Barplot of PINK1 FPKM values faceted by isoforms
}

