\name{TransMeta-Rare-package}
\alias{TransMeta-Rare-package}
\alias{TransMeta-Rare}
\docType{package}
\title{
Trans-ethnic Meta-analysis of Rare Variants in Sequencing Association Studies
}
\description{
This package implements the score test for rare variant associations in trans-ethnic meta-analysis.
}
\details{
\tabular{ll}{
Package: \tab TransMeta-Rare\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2018-06-05\cr
License: \tab GPL (>= 2)\cr
}
This package implements a score test for rare variant associations in trans-ethnic meta-analysis. The method uses summary level score statistics to carry out gene-based meta-analysis for rare variants. 

This package depends on the GMMAT package to implement the logistic mixed model for binary traits in the family-based study design. The GMMAT package is currently only available on Linux platform. It can be installed from GitHub:

library(devtools)

install_github("hanchenphd/GMMAT") 
}
\author{
Jingchunzi Shi, Seunggeun Lee

Maintainer: Jingchunzi (Jing) Shi <shijingc@umich.edu>
}

