\name{TransMeta_Rare_wZ}
\alias{TransMeta_Rare_wZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trans-ethnic meta-analysis with individual level genotype data
}
\description{
A score test for rare variant associations in trans-ethnic meta-analysis. The method uses summary level score statistics to carry out gene-based meta-analysis for rare variants.
}
\usage{
TransMeta_Rare_wZ(Z, obj, combined.weight = FALSE, weights.beta = c(1, 25), 
	Kernel = Kernel, n.Resampling.Copula = 500, Group_Idx = NULL, 
	 rho1 = c(0, 1), rho2 = c(0, 0.09, 0.25, 1), impute.method = "fixed", 
	impute.estimate.maf = 1, missing_cutoff = 0.15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
a numeric genotype matrix with each row as a different individual and each column as a separate snp. Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. Missing genotypes will be imputed using observed MAFs.
}
  \item{obj}{
an output object from the TransMeta_Rare_Null_Model or TransMeta_Rare_Null_Model_EmmaX function.
}
  \item{combined.weight}{
a logical value (default=FALSE) for the type of weighting. If it is TRUE, a weight for each SNP is computed using MAFs that are common across studies. If it is FALSE, group specific weights will be used based on group specific MAFs. 
}
  \item{weights.beta}{
a numeric vector of parameters of beta weights (default=c(1,25)) 
}
  \item{Kernel}{
a kernel matrix which models the correlation structures between the genetic regression coefficients for a one variant across multiple ancestry groups. This kernel matrix can be estimated using the "Get_Kernel_Matrix" function.    %~
}
  \item{n.Resampling.Copula}{
number of iterations for the resampling-based copula algorithm, default value is 500.       %~
}
  \item{Group_Idx}{
a vector of group indicator (default=NULL). If a vector of integers are specified, it assumes causal variants are the same for studies with the same group index, and different for studies with different group indexes. When NULL, studies are assumed to be in different groups with different group indexes.  %~
}
  \item{rho1}{
value of the \eqn{\rho_m}  parameter,  default= c( 0, 1), can also take value 0 or 1. If rho1 =0, it assumes that the population mean genetic effects among the multiple variants are the same; if rho1 = 1, it assumes that the population mean genetic effects among the multiple variants are independently distributed; rho1 = c(0,1) means we do not make prior assumptions regarding the distribution of the population means. %~
}
  \item{rho2}{
value of the \eqn{\rho_k}  parameter,  default= c( 0, 0.09, 0.25, 1).
%% ~
}
  \item{impute.method}{
a method to impute missing genotypes (default= "fixed"). "bestguess" imputes missing genotypes as the most likely values(0,1,2), "random" imputes missing genotypes by generating binomial(2,p) random variables (p = MAF), and "fixed" imputes missing genotypes by assigning the mean genotype value (2p).  %~
}
  \item{impute.estimate.maf}{
a numeric value indicating how to estimate MAFs for the imputation. If impute.estimate.maf=1 (default), MetaSKAT uses study-specific MAFs, in which each study MAFs will be used for the imputation. If impute.estimate.maf=2, all samples in the Z matrix will be used to calculate MAFs for the imputation. Previous versions (< ver 0.6) used impute.estimate.maf=2 as a default.  %~
}
  \item{missing_cutoff}{
a cutoff of the missing rates of SNPs (default=0.15). If the first study has SNPs with missing rates higher than the cutoff, these SNPs in the study will be excluded from the analysis. However, the same SNPs in other studies will not be excluded, if their missing rates are lower than the cutoff. The missing rates are calculated study by study. %~
}
}
\details{
The rows of Z should be matched with phenotypes and covariates. If there are 4 studies, and study 1,2, 3 and 4 have n1, n2, n3 and n4 samples, the first n1, n2, n3, and n4 rows of Z should be genotypes of the first, second, third and forth studies, respectively. Group_Idx is a vector of group index. Suppose the first two studies are European-based, the third study is Asian-based and the last study is African American-based. If you want to run TransMeta-Rare with assuming ancestry group specific heterogeneity, you can set Group_Idx=c(1,1,2,3), which indicates the first two studies belong to the same group. %~
}
\value{
\item{ p.value.normCopula.adj }{p-value for the rare variant association in trans-ethnic meta-analysis}
\item{param }{ estimated parameters of each method.}
%% ...
}

\examples{
############################################################# 

data(y.list)
data(x.list)
data(G.list)
data(K.list)

############################################################# 
# Compute a p-value of the first gene in a population-based study design
Genetic_Kernel = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

obj = TransMeta_Rare_Null_Model(y.list, x.list, n.cohort=4, out_type="C") 
TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj

### groups

Genetic_Kernel_Group = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4, 
	Group_Idx=c(1,1,2,3))

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel_Group, 
	Group_Idx=c(1,1,2,3))$p.value.normCopula.adj

############################################################# 
# Compute a p-value of the first gene in a family-based study design
obj<-TransMeta_Rare_Null_Model_EmmaX(y.list, x.list, n.cohort=4, out_type="C", K.list)

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

