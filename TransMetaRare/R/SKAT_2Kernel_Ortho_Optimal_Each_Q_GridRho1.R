SKAT_2Kernel_Ortho_Optimal_Each_Q_GridRho1 <-
function( Q.all, rho2, rho1, z1.rho1is0, z2.rho1is0, Z1.rho1is0, Z2.rho1is0, z1.rho1is1, z2.rho1is1, Z1.rho1is1, Z2.rho1is1, Phi.tld, Phi.tld.hf,  n.Resampling.Copula){

	n.r<-length(rho2)
	n.q<-dim(Q.all)[1] - n.Resampling.Copula
	n.total = dim(Q.all)[1]

	pval.davis <- matrix(rep(0,2*n.r*n.q),ncol=2*n.r)  ## 2 is the length of rho1 grid, since rho1 = c(0,1)

	for(i in 1:n.r){
		Q.rho1is0<-Q.all[1,i];  Q.rho1is1<-Q.all[1,(i+n.r)]
		if( n.q > 1) Q.boot.rho1is0<-Q.all[2:n.q,i];  Q.boot.rho1is1<-Q.all[2:n.q,(i+n.r)]
		r.corr<-rho2[i]
			
		if(r.corr == 0){ 
			Phi.rho1is0 = t(z1.rho1is0) %*% Phi.tld %*% z1.rho1is0
			Phi.rho1is1 = t(z1.rho1is1) %*% Phi.tld %*% z1.rho1is1
			if( n.q > 1){
				pval.out.rho1is1 = SKAT:::Get_Davies_PVal(Q.rho1is1/2, Phi.rho1is1, Q.boot.rho1is1/2)
				pval.davis[,(i+n.r)] = c(pval.out.rho1is1$p.value, pval.out.rho1is1$p.value.resampling)
			}else{
				pval.davis[,(i+n.r)] = SKAT:::Get_Davies_PVal(Q.rho1is1/2, Phi.rho1is1)$p.value
			}	
		}else if (r.corr == 1){
			Phi.rho1is0 = t(z2.rho1is0) %*% Phi.tld %*% z2.rho1is0
			Phi.rho1is1 = t(z2.rho1is1) %*% Phi.tld %*% z2.rho1is1
			a<- as.matrix(sum(Phi.rho1is1))
			if( n.q > 1){
				pval.out.rho1is1 = SKAT:::Get_Liu_PVal(Q.rho1is1/2, a, Q.boot.rho1is1/2)
				pval.davis[,(i+n.r)] = c(pval.out.rho1is1$p.value, pval.out.rho1is1$p.value.resampling)
			}else{
				pval.davis[,(i+n.r)] = SKAT:::Get_Liu_PVal(Q.rho1is1/2, a)$p.value
			}
		}else{ 
			R.M.rho1is0 = r.corr * Z2.rho1is0 + (1-r.corr)*Z1.rho1is0
			L.rho1is0 = chol(R.M.rho1is0)
			Phi.rho1is0 = L.rho1is0 %*% (Phi.tld %*% t(L.rho1is0))
			R.M.rho1is1 = r.corr * Z2.rho1is1 + (1-r.corr)*Z1.rho1is1
			Phi.rho1is1 = Phi.tld.hf%*%R.M.rho1is1%*%Phi.tld.hf
			if( n.q > 1){
				pval.out.rho1is1 = SKAT:::Get_Davies_PVal(Q.rho1is1/2, Phi.rho1is1, Q.boot.rho1is1/2)
				pval.davis[,(i+n.r)] = c(pval.out.rho1is1$p.value, pval.out.rho1is1$p.value.resampling)
			}else{
				pval.davis[,(i+n.r)] =	SKAT:::Get_Davies_PVal(Q.rho1is1/2, Phi.rho1is1)$p.value
			}
		}

		if( n.q > 1){
			pval.out.rho1is0 = SKAT:::Get_Davies_PVal(Q.rho1is0/2, Phi.rho1is0, Q.boot.rho1is0/2)	
			pval.davis[,i] = c(pval.out.rho1is0$p.value, pval.out.rho1is0$p.value.resampling)
		}else{
			pval.davis[,i] = SKAT:::Get_Davies_PVal(Q.rho1is0/2, Phi.rho1is0)$p.value
		}
		
	}
	
	pmin<-apply(pval.davis,1,min)

	#Sigma = cor( Q.all[(n.q+1):n.total,] ,method = "kendall");
	Sigma = cor( Q.all[(n.q+1):n.total,] ,method = "spearman");

	out<-list(pmin=pmin,pval=pval.davis, Sigma = Sigma)
	return(out)

}
