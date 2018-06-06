\name{G.list}
\alias{G.list}
\docType{data}
\title{Example dataset, a list of genotype matrices}
\description{
a list object of genotypes of all samples. It has 50 elements for 50 genes. Each element is an nxp matrix with n being the total sample size (12,000) over the 4 cohorts and p being the number of SNPs defined in a gene region.  %~
}
\usage{data(G.list)}
\examples{
data(G.list)

length(G.list)

dim(G.list[[1]])

colMeans(G.list[[1]])/2
}
\keyword{datasets}

