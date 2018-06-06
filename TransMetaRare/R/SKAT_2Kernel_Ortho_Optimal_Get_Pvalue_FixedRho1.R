SKAT_2Kernel_Ortho_Optimal_Get_Pvalue_FixedRho1 <-
function(Q.all, z1.1 , z2.1 , rho2, rho1, z1, z2, Z1, Z2, Phi.tld, Phi.tld.hf, n.Resampling.Copula){

	n.r<-length(rho2)
	n.q<-dim(Q.all)[1]-n.Resampling.Copula
	p.m<-dim(z1.1)[2]

	
	# Get Mixture param 
	Each_Info<-SKAT_2Kernel_Ortho_Optimal_Each_Q_FixedRho1(Q.all, rho2, rho1, z1, z2, Z1, Z2, Phi.tld, Phi.tld.hf, n.Resampling.Copula)
	p.val.each <- Each_Info$pval
	pval.copula <- matrix(0, nrow = n.q, ncol =1)
	
	z1.lamb = Get_Lambda(t(z1.1) %*% z1.1)
	if(z1.lamb$zNotZero){
		for(i in 1:n.q){
			pval.copula[i,]<-SKAT_2Kernel_Copula_Optimal_PValue(p.val.each[i,], Each_Info$Sigma)
		}		
	} else {
		pval.copula = Each_Info$pmin
	}
	
	
	multi<-3
	if(length(rho2) < 3){
		multi<-2
	}

	for(i in 1:n.q){
		pval.each<-Each_Info$pval[i,]
		IDX<-which(pval.each > 0)
		
		pval1<-min(pval.each) * multi
		if(pval.copula[i,1] <= 0 || length(IDX) < length(rho2)){
			pval.copula[i,]<-pval1
		}
		
		# if pval==0, use nonzero min each.pval as p-value
		if(pval.copula[i,1] == 0){
			if(length(IDX) > 0){
				pval.copula[i, ] = min(pval.each[IDX])
			}
		}

	
	}
	

	return(list(p.value.copula=pval.copula, p.val.each=p.val.each, Sigma = Each_Info$Sigma))


}
