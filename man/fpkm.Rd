\name{fpkm-methods}
\docType{methods}
\alias{fpkm}
\alias{fpkm-methods}
\alias{fpkm,CuffData-method}
\alias{fpkm,CuffFeature-method}
\alias{fpkm,CuffFeatureSet-method}
\title{
Retrieve FPKM values
}
\description{
Returns a data.frame from @FPKM slot
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
Returns a data.frame of FPKM values.
}
\value{
A data.frame of FPKM-level values for a set of features.
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
	fpkm(PINK1)
}