\name{expressionPlot}

\alias{expressionPlot}
\alias{expressionPlot,CuffFeature-method}
\alias{expressionPlot,CuffFeatureSet-method}

\title{
Expression Plot
}
\description{
A line plot (optionally with confidence intervals) detailing FPKM expression levels across conditions for a given gene(s) or feature(s)
}
\usage{
\S4method{expressionPlot}{CuffFeature}(object, logMode=FALSE, pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE, showStatus=TRUE, replicates=FALSE, ...)
}

\arguments{
  \item{object}{
An object of class ('CuffFeature' or 'CuffGene')
}
  \item{logMode}{
A logical value to draw y-axis (FPKM) on log-10 scale. Default = FALSE.
}
  \item{pseudocount}{
A numeric value added to FPKM to avoid errors on log-10 transformation.
}
  \item{drawSummary}{
A logical valuie. Draws a 'summary' line with mean FPKM values for each condition.
}
  \item{sumFun}{
Function used to determine values for summary line. Default = mean_cl_boot
}
  \item{showErrorbars}{
A logical value whether or not to draw error bars.
}
  \item{showStatus}{
A logical value whether or not to draw visual queues for quantification status of a given gene:condition. Default = TRUE
}
  \item{replicates}{
A logical value whether or not to plot individual replicates or aggregate condition values.
}
  \item{\dots}{
Additional arguments
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
None
}

\examples{
	data(sampleData)
	PINK1 # sample CuffFeature object
	expressionPlot(PINK1) #Line plot of PINK1 FPKM values
	expressionPlot(PINK1@isoforms) #Line plot of PINK1 FPKM values faceted by isoforms
}

