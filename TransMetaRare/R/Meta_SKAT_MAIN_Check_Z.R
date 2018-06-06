Meta_SKAT_MAIN_Check_Z <-
function(Z, n, id_include, SetID, is_dosage=FALSE, impute.method="fixed"){


	m = ncol(Z)
	# OUT

	#############################################
	# Check parameters

	if (class(Z)!= "matrix") stop("Z is not a matrix")
	if (nrow(Z)!=n) stop("Dimensions of y and Z do not match")
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}

	##############################################
	# Check Missing 

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	
	missing.ratio<-rep(0,m)
	for(i in 1:m){
		missing.ratio[i]<-length(which(is.na(Z[,i])))/n
		
		# if missing rate is too high, treat it as unobserved.
		if(missing.ratio[i] > 0.9){
			Z[,i]<-0
		}
	}

	##################################################################
	# doing imputation

	MAF<-colMeans(Z, na.rm = TRUE)/2
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){

		msg<-sprintf("Genotypes of some variants are not the number of minor alleles!")
		#warning(msg,call.=FALSE)
	}
	Z<-SKAT:::Impute(Z,impute.method)

	###########################################
	# Check missing of y and X
	Z<-cbind(Z)
	#id_include1<<-id_include	
	Z.test<-as.matrix(Z[id_include,])

	return(list(Z.test=Z.test, missing_rate=missing.ratio, MAF=MAF , return=0) )

}
