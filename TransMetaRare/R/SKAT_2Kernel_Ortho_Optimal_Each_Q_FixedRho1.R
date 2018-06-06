SKAT_2Kernel_Ortho_Optimal_Each_Q_FixedRho1 <-
function( Q.all, rho2, rho1, z1, z2, Z1, Z2, Phi.tld, Phi.tld.hf, n.Resampling.Copula){

	n.r<-length(rho2)
	n.q<-dim(Q.all)[1] - n.Resampling.Copula
	n.total = dim(Q.all)[1]

	pval.davis<-matrix(rep(0,n.r*n.q),ncol=n.r)  ## 1 is the length of rho1 grid, since rho1 = 0 or 1

	if(rho1 == 0){
		for(i in 1:n.r){
			Q<-Q.all[1,i]; Q.boot<-Q.all[2:n.q,i]
			r.corr<-rho2[i]
			
			if(r.corr == 0){ 
				Phi = t(z1) %*% Phi.tld %*% z1
			}else if (r.corr == 1){
				Phi = t(z2) %*% Phi.tld %*% z2
			}else{ 
				R.M = r.corr * Z2 + (1-r.corr)*Z1
				L = chol(R.M)
				Phi = L %*% (Phi.tld %*% t(L))
			}
			pval.out = SKAT:::Get_Davies_PVal(Q/2, Phi, Q.boot/2)	
			pval.davis[,i] = c(pval.out$p.value, pval.out$p.value.resampling)	
		}
	}else{
		for(i in 1:n.r){
			Q<-Q.all[1,i]; Q.boot<-Q.all[2:n.q,i]
			r.corr<-rho2[i]
			
			if(r.corr == 0){ 
				Phi = t(z1) %*% Phi.tld %*% z1
				pval.out = SKAT:::Get_Davies_PVal(Q/2, Phi, Q.boot/2)
				pval.davis[,i] = c(pval.out$p.value, pval.out$p.value.resampling)		
			}else if (r.corr == 1){
				Phi = t(z2) %*% Phi.tld %*% z2
				a<- as.matrix(sum(Phi))
		            pval.out = SKAT:::Get_Liu_PVal(Q/2, a, Q.boot/2)
				pval.davis[,i] = c(pval.out$p.value, pval.out$p.value.resampling)
			}else{ 
				R.M = r.corr * Z2 + (1-r.corr)*Z1
				#L = chol(R.M)
				#Phi = L %*% (Phi.tld %*% t(L))
				Phi = Phi.tld.hf%*%R.M%*%Phi.tld.hf
				pval.out = SKAT:::Get_Davies_PVal(Q/2, Phi, Q.boot/2)
				pval.davis[,i] = c(pval.out$p.value, pval.out$p.value.resampling)		
			}	
		}
	}	
	pmin<-apply(pval.davis,1,min)

	#Sigma = cor( Q.all[(n.q+1):n.total,] ,method = "kendall");
	Sigma = cor( Q.all[(n.q+1):n.total,] ,method = "spearman");

	out<-list(pmin=pmin,pval=pval.davis, Sigma = Sigma)
	return(out)

}
