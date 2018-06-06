Beta.Weights <-
function(MAF,weights.beta, Cutoff=1, Is.MAF=TRUE){

	n<-length(MAF)
	weights<-rep(0,n)
	Sign<-rep(1,n)
	#print(MAF)
	
	IDX1<-which(MAF > 0.5)
	if(length(IDX1) > 0){
		Sign[IDX1]<--1
		MAF[IDX1]<-1-MAF[IDX1]
	}
	 

	
	IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
	if(length(IDX_0) == n){
		#stop("No polymorphic SNPs")
		weights<-rep(0,n)
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}

	weights = weights * Sign	
	#if(!Is.MAF){
	#	weights1<<-weights
	#	MAF1<<-MAF
	#} else {
	#	weights2<<-weights
	#	MAF2<<-MAF
	#}
	

	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}
