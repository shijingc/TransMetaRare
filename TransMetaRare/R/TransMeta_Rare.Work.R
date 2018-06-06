TransMeta_Rare.Work <-
function(re, n.g, combined.weight=FALSE, n1=NULL, weights.beta=c(1,25), Kernel = Kernel, 
n.Resampling.Copula = n.Resampling.Copula, rho1 = rho1, rho2 = rho2, Group_Idx=NULL, MAF.cutoff=1, Is.MAF=TRUE, missing_cutoff=0.15){


	set.seed(12345)

	# Combined MAF
	p<-length(re[[1]]$MAF)
	
	MAF.Combine=0
	MAF.Groups<-list()
	Map.Groups<-rep(0,n.g)



	MAF.list<-list()
	for(i in 1:n.g){
	
		MAF.list[[i]]<-re[[i]]$MAF
		
	}

	for(i in 1:n.g){
	
		MAF.Combine = MAF.Combine + MAF.list[[i]] * n1[i] / sum(n1)
	}
	
	# If MAF.Combined==0 for all SNP, return p-value 1
	if(sum(MAF.Combine) == 0){
		warning("No polymorphic SNPs!",call.=FALSE)
		return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
	}
	

	# Get MAF.Groups when Group_Idx != NULL
	ID.Groups = unique(Group_Idx)	
	for(j in 1:length(ID.Groups)){
		MAF.Groups[[j]] = 0;
		temp<-which(Group_Idx == ID.Groups[j])
		Map.Groups[temp]<-j
		for(i in temp){
			MAF.Groups[[j]] = MAF.Groups[[j]] + MAF.list[[i]] * n1[i] / sum(n1[temp])
		}
	}


	for(i in 1:n.g){
		if(combined.weight == TRUE){        ### The weight matrix based on ancestry-specific or common MAF is handled here for each study
			weight1<-Beta.Weights(MAF.Combine,weights.beta, MAF.cutoff, Is.MAF=Is.MAF)
		} else {
			j<-Map.Groups[i]
			weight1<-Beta.Weights(MAF.Groups[[j]],weights.beta, MAF.cutoff, Is.MAF=Is.MAF)
		} 
		
		# if missing < missing_cutoff
		
		idx_missing_exclude<-which(re[[i]]$MissingRate > missing_cutoff)
		if(length(idx_missing_exclude) > 0){
			weight1[idx_missing_exclude]<-0
		}
		
		re[[i]]$Score =  re[[i]]$Score * weight1
		re[[i]]$SMat.Summary =  t(t(re[[i]]$SMat.Summary * weight1) * weight1)

		if(!is.null(re[[i]]$Score.Resampling)){
			re[[i]]$Score.Resampling =  re[[i]]$Score.Resampling * weight1
		}
	}


	re.score<-TransMeta_Rare.Work.Groups(re, n.g, ID.Groups, Group_Idx)

	m = re.score$m;   k = re.score$k
	S.tld = re.score$Score;  S.boot.tld = re.score$Score.Resampling;  Phi.tld = re.score$SMat.Summary;  
	res.out = matrix(rnorm(m*k*n.Resampling.Copula), nrow = m*k, ncol = n.Resampling.Copula)

	J.m = rep(1, m);	J.k = rep(1, k);    I.m = diag(1, m);	I.k = diag(1, k)
	M = J.k %*% solve(t(J.k) %*% J.k) %*% t(J.k)

	if( length(rho1)==1 ){
		if(rho1 == 0){
			Z2 = kronecker(J.k %*% t(J.k), I.m)
			Z1 = kronecker((Kernel - M %*% Kernel) %*% (I.k - M), I.m)

			z2 = kronecker(J.k, I.m)	
			L = t(chol(Kernel, pivot = FALSE))
			z1 = kronecker(L - M %*% L , I.m)		
			re = TransMeta_Rare_Optimal_FixedRho1(S.tld, S.boot.tld, Phi.tld, res.out, n.Resampling.Copula, z1, z2, Z1, Z2, rho1, rho2)
		}else{
			Z2 = kronecker(J.k %*% t(J.k), J.m %*% t(J.m))
			Z1 = kronecker((Kernel - M %*% Kernel) %*% (I.k - M), J.m %*% t(J.m))

			z2 = kronecker(J.k, J.m)	
			L = t(chol(Kernel, pivot = FALSE))
			z1 = kronecker(L - M %*% L , J.m)			
			re = TransMeta_Rare_Optimal_FixedRho1(S.tld, S.boot.tld, Phi.tld, res.out, n.Resampling.Copula, z1, z2, Z1, Z2, rho1, rho2)
		} 
	}else{
		Z2.rho1is0 = kronecker(J.k %*% t(J.k), I.m)
		Z1.rho1is0 = kronecker((Kernel - M %*% Kernel) %*% (I.k - M), I.m)
		Z2.rho1is1 = kronecker(J.k %*% t(J.k), J.m %*% t(J.m))
		Z1.rho1is1 = kronecker((Kernel - M %*% Kernel) %*% (I.k - M), J.m %*% t(J.m))

		L = t(chol(Kernel, pivot = FALSE))
		z2.rho1is0 = kronecker(J.k, I.m)
		z1.rho1is0 = kronecker(L - M %*% L , I.m)
		z2.rho1is1 = kronecker(J.k, J.m)
		z1.rho1is1 = kronecker(L - M %*% L , J.m)
		re = TransMeta_Rare_Optimal_GridRho1(S.tld, S.boot.tld, Phi.tld, res.out, n.Resampling.Copula, z1.rho1is0, z2.rho1is0, z1.rho1is1, z2.rho1is1, Z1.rho1is0, Z2.rho1is0, Z1.rho1is1, Z2.rho1is1, rho1, rho2)
	}

	
	return(re)


	
}
