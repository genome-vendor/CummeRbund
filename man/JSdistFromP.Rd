\name{JSdistFromP}
\alias{JSdistFromP}

\title{
Jensen-Shannon distance on rows from a pre-defined vector of probabilities
}
\description{
JSdist takes a matrix of expression probabilites (calculated directly or output from makeprobs()) and returns a matrix of Jensen-Shannon distances between individual rows
and a specific vector of probabilities (q)}
\usage{
JSdistFromP(mat,q)
}

\arguments{
  \item{mat}{
A matrix of expression probabilities (e.g. from makeprobs())
}
  \item{q}{
A vector of expression probabilities.
}
}
\details{
Returns Jensen-Shannon distance for each row of a matrix of probabilities against a provided probability distribution (q)
}
\value{
A vector of JS distances
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
q<-c(100,4,72,8,19)
q<-q/sum(q)
js.distance<-JSdistFromP(mat,q)

}

\keyword{ Jensen-Shannon }
\keyword{ distance }
\keyword{ metric }
\keyword{ JS }