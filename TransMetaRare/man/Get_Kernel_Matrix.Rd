\name{Get_Kernel_Matrix}
\alias{Get_Kernel_Matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to estimate the kernel matrix}
\description{
Use this function to estimate the kernel matrix which models the correlation structures between the genetic regression coefficients for a one variant across multiple ancestry groups.%~
}
\usage{
Get_Kernel_Matrix(y.list, x.list, G, n.cohort, Group_Idx = NULL, 
	is.Genetic.Similarity = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.list}{
a list object for phenotypes. Each element should be a vector of phenotypes. If you have 4 cohorts, it should have 4 elements.
 % ~~
}
  \item{x.list}{
a list object for covariates. Each element should be a vector or a matrix of covariates. If there are 4 cohorts, it should have 4 elements.  %~
}
  \item{G}{
a list object of genotypes of all samples. The number of elements contained in this list should be the total number of gene regions for testing. For example, if testing 50 genes for their associations with the phenotype, then G should have 50 elements for those 50 genes. Each element is an nxp matrix with n being the total sample size over all the cohorts and p being the number of SNPs defined in a gene region. %~
}
  \item{n.cohort}{
a numeric value of the number of cohort.
%~
}
  \item{Group_Idx}{
a vector of group indicator (default=NULL). If a vector of integers are specified, it assumes causal variants are the same for studies with the same group index, and different for studies with different group indexes. When NULL, studies are assumed to be in different groups with different group indexes.  %~
}
  \item{is.Genetic.Similarity}{
to estimate the genetic similarity kernel or group-wise independence kernel. The default value is TRUE, which tells the function to return the genetic similarity kernel; when the value is FALSE, the function returns the group-wise independence kernel.%~
}
}
\examples{
### Genetic similarity kernel
data(y.list)
data(x.list)
data(G.list)


Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

### Genetic similarity kernel with groups
Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4, Group_Idx=c(1,1,2,3))

### Group-wise independence kernel
Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4, is.Genetic.Similarity = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

