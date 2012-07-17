\name{csVolcano}
\alias{csVolcano}
\alias{csVolcano,CuffData-method}
\alias{csVolcano,CuffFeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Volcano Plot
}
\description{
Creates a volcano plot of log fold change in expression vs -log(pval) for a pair of samples (x,y)
}
\usage{
\S4method{csVolcano}{CuffData}(object, x, y, alpha=0.05, showSignificant=TRUE,features=FALSE, xlimits = c(-20, 20), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class CuffData, CuffFeatureSet, or CuffGeneSet
}
  \item{x}{
Sample name from 'samples' table for comparison
}
  \item{y}{
Sample name from 'samples' table for comparison
}
  \item{alpha}{
Provide an alpha cutoff for visualizing significant genes
}
  \item{showSignificant}{
A logical value whether or not to distinguish between significant features or not (by color).
}
  \item{features}{
Will include all fields from 'features' slot in returned ggplot object.
Useful for further manipulations of plot object using feature-level attributes (e.g. gene_type, class_code, etc)
}
  \item{xlimits}{
Set boundaries for x limits to avoid infinity plotting errors. [Default c(-20,20)]
}
  \item{\dots}{
Additional arguments
}
}
\details{
This creates a 'volcano' plot of fold change vs. significance for a pairwise comparison of genes or features across two different samples.
}
\value{
A ggplot2 plot object
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None
}

\examples{
a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object
genes<-a@genes #Create cuffData object for all genes
v<-csVolcano(genes,"hESC","Fibroblasts") # Volcano plot of all genes for conditions x='hESC' and y='Fibroblast'
v #print plot
}

