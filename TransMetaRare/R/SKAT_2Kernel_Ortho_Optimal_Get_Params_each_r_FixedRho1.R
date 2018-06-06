SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r_FixedRho1 <-
function(z1.1, z2.1, rho2){

	c1<-matrix(rep(0,4* length(rho2)), ncol=length(rho2))
	
	A1<-t(z1.1) %*% z1.1
	B1<-t(z2.1) %*% z2.1
	
	A2<-A1 %*% A1
	B2<-B1 %*% B1
	
	A11<-t(z1.1) %*% z2.1
	A22<-A11 %*% t(A11)
	B22<-t(A11) %*% A11
	B333<-t(A11) %*%  A1 %*% A11
	
	#####################################
	#
	
	c1[1,]<-sum(z1.1^2) * (1-rho2) + sum(z2.1^2) * rho2
	c1[2,]<-sum(A1^2) * (1-rho2)^2 + sum(B1^2) * (rho2)^2 + sum(A11^2) * 2 * (1-rho2) * rho2
	c1[3,]<-sum(A2 * A1) * (1-rho2)^3 + sum(B2 * B1) * (rho2)^3 + sum(A22 * A1) * 3 * (1-rho2)^2 * rho2 + sum(B1 * B22) * 3 * (1-rho2) * rho2^2
	c1[4,]<-sum(A2 * A2) * (1-rho2)^4 + sum(B2 * B2) * (rho2)^4 + sum(A22 * A2) * 4 * (1-rho2)^3 * rho2 + sum(B2 * B22) * 4 * (1-rho2) * rho2^3 + sum(B1 * B333) * 4 * (1-rho2)^2 * rho2^2 + sum(B22 * B22) * 2 * (1-rho2)^2 * rho2^2
				
	return(c1)

}
