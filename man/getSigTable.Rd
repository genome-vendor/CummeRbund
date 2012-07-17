\name{getSigTable}
\alias{getSigTable}
\alias{getSigTable,CuffSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getSigTable
}
\description{
Returns the identifiers of significant genes in  a test table - like format.
}
\usage{
\S4method{getSigTable}{CuffSet}(object,alpha=0.05,level='genes')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A CuffSet object (e.g. cuff)
}
  \item{alpha}{
An alpha value by which to filter multiple-testing corrected q-values to determine significance
}
  \item{level}{
Feature level to be queried for significance (must be one of c('genes','isoforms','TSS','CDS')
}
}
\details{
This is a convenience function to quickly retrieve lists of identifiers for genes or features that were determined to be significantly regulated between conditions by cuffdiff.
This function only returns tracking IDs that correspond to tests with an 'OK' status from cuffdiff, NOTEST values are ignored or reported as NA.
By default getSig returns a table of genes x conditions, where the column names represent the pairwise comparisons from the cuffdiff analysis. The values in the table are 1 for features that are significant 
for this comparison and 0 for genes that are not, any failed tests are reported as <NA>.Only includes the features that are significant in at least one comparison.
}		
\value{
A data.frame of pairwise test results.
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
	mySigTable<-getSigTable(a,alpha=0.05,level='genes')
	head(mySigTable)
	
}

