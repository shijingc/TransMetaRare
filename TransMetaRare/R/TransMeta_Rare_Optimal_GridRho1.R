TransMeta_Rare_Optimal_GridRho1 <-
function(S.tld, S.boot.tld = NULL, Phi.tld,  res.out = NULL, n.Resampling.Copula = 500, z1.rho1is0, z2.rho1is0, z1.rho1is1, z2.rho1is1, Z1.rho1is0, Z2.rho1is0, Z1.rho1is1, Z2.rho1is1, rho1, rho2){
		
	Phi.tld.hf = matrix.power(Phi.tld, 0.5)
	z1.rho1is0.1 = Phi.tld.hf %*% z1.rho1is0;	    z2.rho1is0.1 = Phi.tld.hf %*% z2.rho1is0
	z1.rho1is1.1 = Phi.tld.hf %*% z1.rho1is1;	    z2.rho1is1.1 = Phi.tld.hf %*% z2.rho1is1

	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-SKAT_2Kernel_Optimal_Get_Q_GridRho1(S.tld, S.boot.tld, z1.rho1is0, z2.rho1is0, z1.rho1is1, z2.rho1is1, z1.rho1is0.1, z2.rho1is0.1, z1.rho1is1.1, z2.rho1is1.1, rho2, n.Resampling.Copula, res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.boot, out.Q$Q.r.res)  


	##################################################
	# Compute P-values 
	#################################################

	out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue_GridRho1(Q.all, z1.rho1is0.1, z2.rho1is0.1, z1.rho1is1.1, z2.rho1is1.1,  z1.rho1is0, z2.rho1is0, Z1.rho1is0, Z2.rho1is0,  z1.rho1is1, z2.rho1is1, Z1.rho1is1, Z2.rho1is1, Phi.tld, Phi.tld.hf, rho1, rho2,  n.Resampling.Copula)

	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each   ## p-value of U.tau for each rho2 value
	n.q<-dim(Q.all)[1] - n.Resampling.Copula
	param$Q.each<-Q.all[1:n.q,]  ## value of U.tau for each rho2 value 
	param$rho<-cbind(rep(c(0,1), each = 4), rep(rho2, 2))          ## grid value of rho1, rho2 combination
	param$minp<- apply(param$p.val.each, 1, 'min')   ## value of test stat. T  (T = min p-value of U.tau for each rho2 value)
	param$Sigma = out$Sigma

	id_temp<- apply(as.matrix(seq(from = 1, to = n.q)), 1, function(x) which(param$p.val.each[x,] == param$minp[x])[1] )
	param$rho_est<-param$rho[id_temp,]   ## the (rho1, rho2) value whose corresponding U.tau value yields the smallest p-value (i.e. T equals that p-value)

	#p.value <- out$p.value[1]  ## p-value of the test stat. T
	p.value.normCopula = out$p.value.copula[,1]
	p.value.normCopula.adj <-  apply( as.matrix(seq(from = 1, to = n.q)), 1, function(x) min(out$p.value.copula[x,1], max(out$p.val.each[x,])) )


  	re<-list(p.value.normCopula = p.value.normCopula, p.value.normCopula.adj = p.value.normCopula.adj, param=param )  
  	
	return(re)	

}
