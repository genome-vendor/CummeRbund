\name{repFpkm-methods}
\docType{methods}
\alias{repFpkm}
\alias{repFpkm-methods}
\alias{repFpkm,CuffData-method}
\alias{repFpkm,CuffFeature-method}
\alias{repFpkm,CuffFeatureSet-method}
\title{
Retrieve FPKM values
}
\description{
Returns a data.frame from @repFpkm slot
}
\section{Methods}{
\describe{

\item{\code{signature(object = "CuffData")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "CuffFeature")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "CuffFeatureSet")}}{
%%  ~~describe this method here~~
}
}}

\details{
Returns a data.frame of replicate FPKM values and associated statistics.
}
\value{
A data.frame of replicate-level FPKM values and associated statistics for a set of features.
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
	repFpkm(PINK1)
}