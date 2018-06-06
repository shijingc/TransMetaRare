\name{x.list}
\alias{x.list}
\docType{data}
\title{Example dataset, a list of non-genetic adjusting covariate matrices}
\description{
a list object of covariates. It has 4 elements for 4 study cohorts. Each element is a matrix of covariates. The first, third and last elements have two covariates (two columns), and the second element has one covariate (one column).
}
\usage{data(x.list)}
\examples{
data(x.list)

length(x.list)

head(x.list[[1]])
}
\keyword{datasets}

