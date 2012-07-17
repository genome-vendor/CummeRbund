\name{getSig}
\alias{getSig}
\alias{getSig,CuffSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getSig
}
\description{
Returns the identifiers of significant genes in a vector format.
}
\usage{
\S4method{getSig}{CuffSet}(object,x,y,alpha=0.05,level='genes')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A CuffSet object (e.g. cuff)
}
  \item{x}{
Optional argument to restrict significance results to one pairwise comparison.  Must be used with a 'y' argument to specificy the other half of the pair.
}
  \item{y}{
See x.
}
  \item{alpha}{
An alpha value by which to filter multiple-testing corrected q-values to determine significance
}
  \item{level}{
Feature level to be queried for significance (must be one of c('genes','isoforms','TSS','CDS')
}

}
\details{
This is a convenience function to quickly retrieve vectors of identifiers for genes or features that were determined to be significantly regulated between conditions by cuffdiff.
This function only returns tracking IDs that correspond to tests with an 'OK' status from cuffdiff, NOTEST values are ignored.
By default getSig returns a vector of IDs for all pairwise comparisons together. If you specify an 'x' AND 'y' values as sample names, then only the features that are significant in that
particular pairwise comparison are reported, after appropriate multiple testing correction of output p-values.
}		
\value{
A vector of feature IDs.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}


\examples{

	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data in sample directory and creates CuffSet object 'a'
	mySig<-getSig(a,x='hESC',y='Fibroblasts',alpha=0.05,level='genes')
	head(mySig)
	
}

