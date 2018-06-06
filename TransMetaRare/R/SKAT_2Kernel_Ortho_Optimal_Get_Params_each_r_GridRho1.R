SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r_GridRho1 <-
function(z1.rho1is0.1, z2.rho1is0.1, z1.rho1is1.1, z2.rho1is1.1, rho2){

	c1.rho1is0<-matrix(rep(0,4* length(rho2)), ncol=length(rho2))
	
	A1.rho1is0 <-t(z1.rho1is0.1) %*% z1.rho1is0.1
	B1.rho1is0 <-t(z2.rho1is0.1) %*% z2.rho1is0.1
	
	A2.rho1is0 <-A1.rho1is0 %*% A1.rho1is0
	B2.rho1is0 <-B1.rho1is0 %*% B1.rho1is0
	
	A11.rho1is0 <-t(z1.rho1is0.1) %*% z2.rho1is0.1
	A22.rho1is0<-A11.rho1is0 %*% t(A11.rho1is0)
	B22.rho1is0<-t(A11.rho1is0) %*% A11.rho1is0
	B333.rho1is0<-t(A11.rho1is0) %*%  A1.rho1is0 %*% A11.rho1is0
	
	c1.rho1is0[1,]<-sum(z1.rho1is0.1^2) * (1-rho2) + sum(z2.rho1is0.1^2) * rho2
	c1.rho1is0[2,]<-sum(A1.rho1is0^2) * (1-rho2)^2 + sum(B1.rho1is0^2) * (rho2)^2 + sum(A11.rho1is0^2) * 2 * (1-rho2) * rho2
	c1.rho1is0[3,]<-sum(A2.rho1is0 * A1.rho1is0) * (1-rho2)^3 + sum(B2.rho1is0 * B1.rho1is0) * (rho2)^3 + sum(A22.rho1is0 * A1.rho1is0) * 3 * (1-rho2)^2 * rho2 + sum(B1.rho1is0 * B22.rho1is0) * 3 * (1-rho2) * rho2^2
	c1.rho1is0[4,]<-sum(A2.rho1is0 * A2.rho1is0) * (1-rho2)^4 + sum(B2.rho1is0 * B2.rho1is0) * (rho2)^4 + sum(A22.rho1is0 * A2.rho1is0) * 4 * (1-rho2)^3 * rho2 + sum(B2.rho1is0 * B22.rho1is0) * 4 * (1-rho2) * rho2^3 + sum(B1.rho1is0 * B333.rho1is0) * 4 * (1-rho2)^2 * rho2^2 + sum(B22.rho1is0 * B22.rho1is0) * 2 * (1-rho2)^2 * rho2^2
	
	
	#####################################
			
	c1.rho1is1<-matrix(rep(0,4* length(rho2)), ncol=length(rho2))
	
	A1.rho1is1<-t(z1.rho1is1.1) %*% z1.rho1is1.1
	B1.rho1is1<-t(z2.rho1is1.1) %*% z2.rho1is1.1
	
	A2.rho1is1<-A1.rho1is1 %*% A1.rho1is1
	B2.rho1is1<-B1.rho1is1 %*% B1.rho1is1
	
	A11.rho1is1<-t(z1.rho1is1.1) %*% z2.rho1is1.1
	A22.rho1is1<-A11.rho1is1 %*% t(A11.rho1is1)
	B22.rho1is1<-t(A11.rho1is1) %*% A11.rho1is1
	B333.rho1is1<-t(A11.rho1is1) %*%  A1.rho1is1 %*% A11.rho1is1
	
	c1.rho1is1[1,]<-sum(z1.rho1is1.1^2) * (1-rho2) + sum(z2.rho1is1.1^2) * rho2
	c1.rho1is1[2,]<-sum(A1.rho1is1^2) * (1-rho2)^2 + sum(B1.rho1is1^2) * (rho2)^2 + sum(A11.rho1is1^2) * 2 * (1-rho2) * rho2
	c1.rho1is1[3,]<-sum(A2.rho1is1 * A1.rho1is1) * (1-rho2)^3 + sum(B2.rho1is1 * B1.rho1is1) * (rho2)^3 + sum(A22.rho1is1 * A1.rho1is1) * 3 * (1-rho2)^2 * rho2 + sum(B1.rho1is1 * B22.rho1is1) * 3 * (1-rho2) * rho2^2
	c1.rho1is1[4,]<-sum(A2.rho1is1 * A2.rho1is1) * (1-rho2)^4 + sum(B2.rho1is1 * B2.rho1is1) * (rho2)^4 + sum(A22.rho1is1 * A2.rho1is1) * 4 * (1-rho2)^3 * rho2 + sum(B2.rho1is1 * B22.rho1is1) * 4 * (1-rho2) * rho2^3 + sum(B1.rho1is1 * B333.rho1is1) * 4 * (1-rho2)^2 * rho2^2 + sum(B22.rho1is1 * B22.rho1is1) * 2 * (1-rho2)^2 * rho2^2
		
	c1 = list();	c1$c1.rho1is0 = c1.rho1is0;	c1$c1.rho1is1 = c1.rho1is1
	return(c1)

}
