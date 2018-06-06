MAC <-
function(MAF.list, n.ancestry){ ##given a gene, count number of shared variants btwn two ancestry groups, and total number of varints in the gene
	SNP = list()
	for(l in 1:n.ancestry){
		SNP[[l]] = which(MAF.list[,l] > 0);  
	}

	MAC.num =  matrix(0, nrow = n.ancestry, ncol = n.ancestry)
	MAC.den =  matrix(dim(MAF.list)[1], nrow = n.ancestry, ncol = n.ancestry)

	for(i in 1:(n.ancestry-1)){
		for(j in (i+1):n.ancestry){
			if( length(union( SNP[[i]], SNP[[j]])) > 0 ){
				MAC.num[i,j] = MAC.num[j,i] = length(intersect( SNP[[i]], SNP[[j]]))
			}else{
				MAC.num[i,j] = MAC.num[j,i] = 0
			}
		}
	}
	
	rs = list()
	rs$MAC.num = MAC.num;  rs$MAC.den = MAC.den
	return(rs)
}
