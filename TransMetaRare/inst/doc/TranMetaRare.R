### R code from vignette source 'TranMetaRare.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: data
###################################################
###########################################################################
####### To run this code, first download and unzip example files
###########################################################################

library(TransMetaRare)
y1 = as.matrix(read.table('y1.txt'))
y2 = as.matrix(read.table('y2.txt'))
y3 = as.matrix(read.table('y3.txt'))
y4 = as.matrix(read.table('y4.txt'))
y.list = list() 	
y.list[[1]] = y1 	
y.list[[2]] = y2	
y.list[[3]] = y3	
y.list[[4]] = y4;
#########################

x1 = as.matrix(read.table('x1.txt'))
x2 = as.matrix(read.table('x2.txt'))
x3 = as.matrix(read.table('x3.txt'))
x4 = as.matrix(read.table('x4.txt'))
x.list = list() 	
x.list[[1]] = x1 	
x.list[[2]] = as.matrix(x2[,1])
x.list[[3]] = x3	
x.list[[4]] = x4

#########################
Kins_EUR = as.matrix(read.table('Kins_EUR.txt'))
Kins_OW = as.matrix(read.table('Kins_OW.txt'))
K.list = list()  
K.list[[1]] = Kins_EUR
K.list[[2]] = Kins_OW
K.list[[3]] = Kins_OW
K.list[[4]] = Kins_OW

##########################
G.list = list()
for(i in 1:50){
	Z = as.matrix(read.table( paste("Z",i,".txt",sep=""))  )
	G.list[[i]] = Z
}



###################################################
### code chunk number 2: Test1
###################################################
# Compute a p-value of the first gene in a population-based study design

Genetic_Kernel = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4)

obj = TransMeta_Rare_Null_Model(y.list, x.list, n.cohort=4, out_type="C") 

rs1 = TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)

rs1$p.value.normCopula.adj



###################################################
### code chunk number 3: Test2
###################################################
### with groups

Genetic_Kernel_G = Get_Kernel_Matrix(y.list, x.list, G.list, n.cohort = 4,
		 Group_Idx=c(1,1,2,3))

rs2 = TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel_G, 
		Group_Idx=c(1,1,2,3))

rs2$p.value.normCopula.adj



###################################################
### code chunk number 4: Test3
###################################################
# Compute a p-value of the first gene in a family-based study design

obj<-TransMeta_Rare_Null_Model_EmmaX(y.list, x.list, n.cohort=4, out_type="C", K.list)

rs3 = TransMeta_Rare_wZ (G.list[[1]], obj, Kernel = Genetic_Kernel)

rs3$p.value.normCopula.adj






