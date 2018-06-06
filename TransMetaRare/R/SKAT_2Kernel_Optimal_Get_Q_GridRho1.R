SKAT_2Kernel_Optimal_Get_Q_GridRho1 <-
function(S.tld, S.boot.tld = NULL, z1.rho1is0, z2.rho1is0, z1.rho1is1, z2.rho1is1, z1.rho1is0.1, z2.rho1is0.1, z1.rho1is1.1, z2.rho1is1.1, rho2, n.Resampling.Copula, res.out){

	n.r<-length(rho2)  
	
	Q.rho1is0.r <- Q.rho1is1.r <- rep(0,n.r)
	Q.r.res<-NULL
	Q.boot<-NULL	
	
	temp1.rho1is0<-t(S.tld) %*% z1.rho1is0
	temp2.rho1is0<-t(S.tld) %*% z2.rho1is0
	temp1.rho1is1<-t(S.tld) %*% z1.rho1is1
	temp2.rho1is1<-t(S.tld) %*% z2.rho1is1

	for(i in 1:n.r){
		r.corr <- rho2[i]
		Q1.rho1is0 <-  (1-r.corr) * rowSums(temp1.rho1is0^2)
		Q2.rho1is0 <- r.corr * rowSums(temp2.rho1is0^2)
		Q.rho1is0.r[i] <- Q1.rho1is0 + Q2.rho1is0
		
		Q1.rho1is1 <-  (1-r.corr) * rowSums(temp1.rho1is1^2)
		Q2.rho1is1 <- r.corr * rowSums(temp2.rho1is1^2)
		Q.rho1is1.r[i] <- Q1.rho1is1 + Q2.rho1is1

	}
	Q.r = c(Q.rho1is0.r, Q.rho1is1.r)

  	if(n.Resampling.Copula > 0){
	
		temp1.rho1is0<-t(res.out) %*% z1.rho1is0.1
		temp2.rho1is0<-t(res.out) %*% z2.rho1is0.1
		temp1.rho1is1<-t(res.out) %*% z1.rho1is1.1
		temp2.rho1is1<-t(res.out) %*% z2.rho1is1.1

		Q.rho1is0.r.res <- Q.rho1is1.r.res<- matrix(rep(0,n.Resampling.Copula *n.r),ncol=n.r)

		for(i in 1:n.r){
			r.corr<-rho2[i]
			Q1.rho1is0 <-  (1-r.corr) * rowSums(temp1.rho1is0^2)
			Q2.rho1is0 <- r.corr * rowSums(temp2.rho1is0^2)
			Q.rho1is0.r.res[,i] <- Q1.rho1is0 + Q2.rho1is0
		
			Q1.rho1is1 <-  (1-r.corr) * rowSums(temp1.rho1is1^2)
			Q2.rho1is1 <- r.corr * rowSums(temp2.rho1is1^2)
			Q.rho1is1.r.res[,i] <- Q1.rho1is1 + Q2.rho1is1
		}
		Q.r.res = cbind(Q.rho1is0.r.res, Q.rho1is1.r.res)
  	}else{
		Q.r.res = NULL
	}
	

	if(!is.null(S.boot.tld)){

		temp1.rho1is0<-t(S.boot.tld) %*% z1.rho1is0
		temp2.rho1is0<-t(S.boot.tld) %*% z2.rho1is0
		temp1.rho1is1<-t(S.boot.tld) %*% z1.rho1is1
		temp2.rho1is1<-t(S.boot.tld) %*% z2.rho1is1
		n.moments<-dim(S.boot.tld)[2]
		Q.rho1is0.r.boot <- Q.rho1is1.r.boot <- matrix(rep(0,n.moments *n.r),ncol=n.r)

		for(i in 1:n.r){
			r.corr<-rho2[i]
			Q1.rho1is0 <-  (1-r.corr) * rowSums(temp1.rho1is0^2)
			Q2.rho1is0 <- r.corr * rowSums(temp2.rho1is0^2)
			Q.rho1is0.r.boot[,i] <- Q1.rho1is0 + Q2.rho1is0
		
			Q1.rho1is1 <-  (1-r.corr) * rowSums(temp1.rho1is1^2)
			Q2.rho1is1 <- r.corr * rowSums(temp2.rho1is1^2)
			Q.rho1is1.r.boot[,i] <- Q1.rho1is1 + Q2.rho1is1
		}
		
		Q.boot = cbind(Q.rho1is0.r.boot, Q.rho1is1.r.boot)
	}else{
		Q.boot = NULL
	}


	re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.boot=Q.boot)
	return(re)

}
