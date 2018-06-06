Get_Lambda <-
function(K){

	out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)
	#print(out.s$values)

	#out.s1<-eigen(K,symmetric=TRUE)
	#print(out.s1$values)
	
	lambda1<-out.s$values
	IDX1<-which(lambda1 >= 0)

	# eigenvalue bigger than sum(eigenvalues)/1000
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
	#cat("Lambda:", lambda1, "\n")
	#K1<<-K
	
	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	}
	lambda<-lambda1[IDX2]
	
	lambda0 = rep(0, length(lambda))
	IDX3 = which(lambda < 10^(-20)) ## very very small eigenvalues, set those to zero

	if(sum(lambda < 10^(-20)) == length(lambda)){ 
		zNotZero = 0  ## all eigenvalues are zeros (small enough to be view as zero)
	} else {
		zNotZero = 1 ## at least one eigenvaleu is bigger than zero
	}
	
	lambda[IDX3] = lambda0[IDX3] ## set those small eigenvalues to zero

	res = list()
	res$lambda = lambda; res$zNotZero = zNotZero

	return(res)

}
