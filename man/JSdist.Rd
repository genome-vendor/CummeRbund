\name{JSdist}
\alias{JSdist}

\title{
Jensen-Shannon distance on columns
}
\description{
JSdist takes a matrix of expression probabilites (calculated directly or output from makeprobs()) and returns a dist object 
of the pairwise Jensen-Shannon distances between columns
}
\usage{
JSdist(mat)
}

\arguments{
  \item{mat}{
A matrix of expression probabilities (e.g. from makeprobs())
}
}
\details{
Returns pairwise Jensen-Shannon distance (in the form of a dist object) for a matrix of probabilities (by column)
}
\value{
A dist object of pairwise J-S distances between columns.
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

mat<-matrix(sample(1:50,50),10)
probs<-makeprobs(mat)
js.distance<-JSdist(probs)

}

\keyword{ Jensen-Shannon }
\keyword{ distance }
\keyword{ metric }
\keyword{ JS }