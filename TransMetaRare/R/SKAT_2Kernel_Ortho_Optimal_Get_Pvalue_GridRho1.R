SKAT_2Kernel_Ortho_Optimal_Get_Pvalue_GridRho1 <-
function(Q.all, z1.rho1is0.1, z2.rho1is0.1, z1.rho1is1.1, z2.rho1is1.1,  z1.rho1is0, z2.rho1is0, Z1.rho1is0, Z2.rho1is0,  z1.rho1is1, z2.rho1is1, Z1.rho1is1, Z2.rho1is1, Phi.tld, Phi.tld.hf, rho1, rho2,  n.Resampling.Copula){

	n.r<-length(rho2)
	n.q<-dim(Q.all)[1]-n.Resampling.Copula
	p.rho1is0.m<-dim(z1.rho1is0.1)[2];	p.rho1is1.m<-dim(z1.rho1is1.1)[2]

	
	# Get Mixture param 
	Each_Info<-SKAT_2Kernel_Ortho_Optimal_Each_Q_GridRho1(Q.all, rho2, rho1, z1.rho1is0, z2.rho1is0, Z1.rho1is0, Z2.rho1is0, z1.rho1is1, z2.rho1is1, Z1.rho1is1, Z2.rho1is1, Phi.tld, Phi.tld.hf,  n.Resampling.Copula)
	p.val.each <- Each_Info$pval
	pval.copula <- matrix(0, nrow = n.q, ncol =1)
	
	z1.rho1is0.lamb = Get_Lambda(t(z1.rho1is0.1) %*% z1.rho1is0.1)
	z1.rho1is1.lamb = Get_Lambda(t(z1.rho1is1.1) %*% z1.rho1is1.1)

	if(z1.rho1is0.lamb$zNotZero || z1.rho1is1.lamb$zNotZero){
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
