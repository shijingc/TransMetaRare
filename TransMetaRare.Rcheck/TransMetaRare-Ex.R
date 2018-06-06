pkgname <- "TransMetaRare"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('TransMetaRare')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("G.list")
### * G.list

flush(stderr()); flush(stdout())

### Name: G.list
### Title: Example dataset, a list of genotype matrices
### Aliases: G.list
### Keywords: datasets

### ** Examples

data(G.list)

length(G.list)

dim(G.list[[1]])

colMeans(G.list[[1]])/2



cleanEx()
nameEx("Get_Kernel_Matrix")
### * Get_Kernel_Matrix

flush(stderr()); flush(stdout())

### Name: Get_Kernel_Matrix
### Title: A function to estimate the kernel matrix
### Aliases: Get_Kernel_Matrix
### Keywords: ~kwd1 ~kwd2

### ** Examples

### Genetic similarity kernel
data(y.list)
data(x.list)
data(G.list)


Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

### Genetic similarity kernel with groups
Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4, Group_Idx=c(1,1,2,3))

### Group-wise independence kernel
Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4, is.Genetic.Similarity = FALSE)




cleanEx()
nameEx("K.list")
### * K.list

flush(stderr()); flush(stdout())

### Name: K.list
### Title: Example dataset, a list of kinship matrices
### Aliases: K.list
### Keywords: datasets

### ** Examples

data(K.list)

length(K.list)

head(K.list[[1]])



cleanEx()
nameEx("TransMeta_Rare_Null_Model")
### * TransMeta_Rare_Null_Model

flush(stderr()); flush(stdout())

### Name: TransMeta_Rare_Null_Model
### Title: Get parameters and residuals from H0 for population-based study
###   design
### Aliases: TransMeta_Rare_Null_Model

### ** Examples

############################################################# 
data(y.list)
data(x.list)
data(G.list)

# Compute a p-value of the first gene 
Genetic_Kernel = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

obj = TransMeta_Rare_Null_Model(y.list, x.list, n.cohort=4, out_type="C") 

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj

############################################################# 
# If you want to use the intercept-only model for the 2nd cohort 
x.list[[2]]<-"intercept" 

obj = TransMeta_Rare_Null_Model (y.list, x.list, n.cohort=4, out_type="C") 

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj





cleanEx()
nameEx("TransMeta_Rare_Null_Model_EmmaX")
### * TransMeta_Rare_Null_Model_EmmaX

flush(stderr()); flush(stdout())

### Name: TransMeta_Rare_Null_Model_EmmaX
### Title: Get parameters and residuals from H0 for family-based study
###   design
### Aliases: TransMeta_Rare_Null_Model_EmmaX

### ** Examples

############################################################# 

data(y.list)
data(x.list)
data(G.list)
data(K.list)

# Compute a p-value of the first gene 
Genetic_Kernel = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

obj = TransMeta_Rare_Null_Model_EmmaX(y.list, x.list, n.cohort=4, out_type="C", K.list)

TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)$p.value.normCopula.adj





cleanEx()
nameEx("TransMeta_Rare_wZ")
### * TransMeta_Rare_wZ

flush(stderr()); flush(stdout())

### Name: TransMeta_Rare_wZ
### Title: Trans-ethnic meta-analysis with individual level genotype data
### Aliases: TransMeta_Rare_wZ
### Keywords: ~kwd1 ~kwd2

### ** Examples

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






cleanEx()
nameEx("x.list")
### * x.list

flush(stderr()); flush(stdout())

### Name: x.list
### Title: Example dataset, a list of non-genetic adjusting covariate
###   matrices
### Aliases: x.list
### Keywords: datasets

### ** Examples

data(x.list)

length(x.list)

head(x.list[[1]])



cleanEx()
nameEx("y.list")
### * y.list

flush(stderr()); flush(stdout())

### Name: y.list
### Title: Example dataset, a list of phenotype vectors
### Aliases: y.list
### Keywords: datasets

### ** Examples

data(y.list)

length(y.list)

length(y.list[[1]])



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
