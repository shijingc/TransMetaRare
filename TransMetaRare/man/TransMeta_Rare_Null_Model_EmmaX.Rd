\name{TransMeta_Rare_Null_Model_EmmaX}
\alias{TransMeta_Rare_Null_Model_EmmaX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get parameters and residuals from H0 for family-based study design
}
\description{
Compute model parameters and residuals under the null model (H0) of no associations for family-based study design. It can be used only when individual level data are available.
}
\usage{
TransMeta_Rare_Null_Model_EmmaX(y.list, x.list, n.cohort, out_type = "C", 
	K.list = NULL, Kin.File.list = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y.list}{
a list object for phenotypes. Each element should be a vector of phenotypes. If you have 4 cohorts, it should have 4 elements.
}
  \item{x.list}{
a list object for covariates. Each element should be a vector or a matrix of covariates. If there are 4 cohorts, it should have 4 elements. If there are no covariates to adjust for, the element should be "intercept". See the examples.
}
  \item{n.cohort}{
a numeric value of the number of cohort.
}
  \item{out_type}{
an indicator for the outcome type. "C" for continuous outcomes and "D" for dichotomous outcomes. 
}
  \item{K.list}{
a list object of the kinship matrices. If K.list=NULL, the function reads files in Kin.File.list to obtain kinship matrices.
}
  \item{Kin.File.list}{
a list object of emmax-kin output file names.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns an object that has model parameters and residuals. The returned object will be used to run TransMeta_Rare_wZ.
}
\examples{
############################################################# 

data(y.list)
data(x.list)
data(G.list)
data(K.list)

# Compute a p-value of the first gene 
Genetic_Kernel = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

obj = TransMeta_Rare_Null_Model_EmmaX(y.list, x.list, n.cohort=4, out_type="C", K.list)

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj


}

