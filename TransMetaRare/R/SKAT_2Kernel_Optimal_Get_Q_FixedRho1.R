SKAT_2Kernel_Optimal_Get_Q_FixedRho1 <-
function(S.tld, S.boot.tld = NULL, z1, z2, z1.1, z2.1, rho2, n.Resampling.Copula, res.out){

	n.r<-length(rho2)
	
	Q.r<-rep(0,n.r)
	Q.r.res<-NULL
	Q.boot<-NULL	
	
	temp1<-t(S.tld) %*% z1
	temp2<-t(S.tld) %*% z2
	for(i in 1:n.r){
		r.corr<-rho2[i]
		Q1<-  (1-r.corr) * rowSums(temp1^2)
		Q2<- r.corr * rowSums(temp2^2)
		Q.r[i]<-Q1 + Q2
	}

  	if(n.Resampling.Copula > 0){
	
		temp1<-t(res.out) %*% z1.1
		temp2<-t(res.out) %*% z2.1
		Q.r.res<-matrix(rep(0,n.Resampling.Copula *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-rho2[i]
			Q1<-  (1-r.corr) * rowSums(temp1^2)
			Q2<- r.corr * rowSums(temp2^2)
			Q.r.res[,i]<-Q1 + Q2
		}
  	}

	if(!is.null(S.boot.tld)){

		temp1<-t(S.boot.tld) %*% z1
		temp2<-t(S.boot.tld) %*% z2
		n.moments<-dim(S.boot.tld)[2]
		Q.boot<-matrix(rep(0,n.moments *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-rho2[i]
			Q1<- (1-r.corr) * rowSums(temp1^2)
			Q2<- r.corr * rowSums(temp2^2)
			Q.boot[,i]<-Q1 + Q2
		}

	}

	re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.boot=Q.boot)
	return(re)

}
