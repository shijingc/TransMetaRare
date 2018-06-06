est.MAC <-
function(n.gene, n.ancestry, MAF.list.Groups){

	MAC.num = MAC.den = matrix(0, nrow = n.ancestry, ncol = n.ancestry)

	for(g in 1:n.gene){ 

		MAC.temp = MAC(MAF.list.Groups[[g]], n.ancestry) 
		MAC.num = MAC.num + MAC.temp$MAC.num;  MAC.den = MAC.den + MAC.temp$MAC.den;
	}

	MAC = MAC.num/MAC.den
	diag(MAC) = rep(1, n.ancestry)
	return(MAC)
}
