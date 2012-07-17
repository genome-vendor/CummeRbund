\name{count-methods}
\docType{methods}
\alias{count}
\alias{count-methods}
\alias{count,CuffData-method}
\alias{count,CuffFeature-method}
\alias{count,CuffFeatureSet-method}
\title{
Retrieve count values (raw and normalized)
}
\description{
Returns a data.frame from @count slot
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
Returns a data.frame of count values.
}
\value{
A data.frame of count-level values for a set of features.
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
	count(PINK1)
}