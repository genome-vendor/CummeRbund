\name{addFeatures}

\alias{addFeatures}
\alias{addFeatures,CuffSet-method}
\alias{addFeatures,CuffData-method}

\title{
addFeatures
}
\description{
Adds a data.frame of features to a the SQLite backend database.
}
\usage{
\S4method{addFeatures}{CuffSet}(object, features, level="genes", ...)

}

\arguments{
  \item{object}{
An object of class ('CuffSet' or 'CuffData')
}
  \item{features}{
A data.frame of features to add.  1st column MUST contain ids (ie. gene_id for 'gene' features, isoform_id for 'isoform' features, etc)
}

  \item{level}{
One of c('genes','isoforms','TSS','CDS') to indicate which type of featurs you are being added, and to what data-level.
  }

  \item{\dots}{
Additional arguments.
}
}
\details{
None
}
\value{
None
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
	#None yet.
}
