TransMeta_Rare_wZ <-
function(Z, obj, combined.weight=FALSE, weights.beta=c(1,25), Kernel = Kernel, 
 n.Resampling.Copula = 500, rho1 = c(0,1), rho2 = c(0, 0.09, 0.25, 1),  Group_Idx=NULL, impute.method="fixed", impute.estimate.maf=1, missing_cutoff=0.15){


	if(is.matrix(Z)!= TRUE){
		stop("ERROR: Z is not a matrix!")
	}
	if(class(obj)!= "META_NULL_Model" && class(obj)!= "META_NULL_Model_EmmaX"){
		stop("ERROR: obj class is not either META_NULL_Model or META_NULL_Model_EmmaX!")
	}
	
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 
	
	nSNP<-ncol(Z)
	n<-nrow(Z)
	n.g<-obj$n.g
	
	# Calculate missing rate : Since Z is imputed before calling Meta_SKAT_SaveData when impute.estimate.maf==1, 
	# missing rate should be calculated before calling this function.
	
	missing.ratio.Z = matrix(rep(0, nSNP*n.g), nrow=n.g)
	for(i in 1:n.g){
		ID<-obj$ID[[i]]
		Z1<-as.matrix(Z[ID,])
		n1<-nrow(Z1)
		for(j in 1:nSNP){
			missing.ratio.Z[i,j]<-length(which(is.na(Z1[,j])))/n1
		}
	}
	
	Is.impute.cohortwise= FALSE
	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )

		warning(msg,call.=FALSE)
		
		if(impute.estimate.maf==1){
			Is.impute.cohortwise = TRUE
		
		} else if(impute.estimate.maf==2){
			Z<-SKAT:::Impute(Z,impute.method=impute.method)
		} else {
			stop("ERROR: impute.estimate.mat is wrong! it should be either 1 or 2")
		}
	} 
	
	
	re1<-list()
	for(i in 1:n.g){

		ID<-obj$ID[[i]]
		Z1<-as.matrix(Z[ID,])
		re1[[i]]<-Meta_SKAT_SaveData(Z1, obj$out[[i]], SetID=NULL, impute.method = impute.method)
		
		re1[[i]]$MissingRate = missing.ratio.Z[i,]
	}

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.g
	}

	re = TransMeta_Rare.Work(re1, n.g, combined.weight, n1=obj$n.each, weights.beta=weights.beta, Kernel = Kernel, n.Resampling.Copula = n.Resampling.Copula, rho1=rho1, rho2 = rho2, 
	, Group_Idx=Group_Idx, missing_cutoff=missing_cutoff)

	return(re)

}
