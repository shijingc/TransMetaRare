Get_Kernel_Matrix <-
function(y.list, x.list, G, n.cohort, Group_Idx=NULL, is.Genetic.Similarity = TRUE ){

	n.g = n.cohort
	n.gene = length(G)  ## total number of gene regions for testing

	if(is.null(Group_Idx)){
		Group_Idx<-1:n.g
	}

	ID.position<-list()
	n.each<-rep(0,n.cohort)
	count<-1
	n.x = length(x.list)
	for(i in 1:n.cohort){
		y = y.list[[i]]
		X = x.list[[i]]
		n1<-length(y.list[[i]])
		ID.position[[i]]<-count:(count +n1-1)
		count<-count + n1

		n.each[i]<-n1
	}


	MAF.list<-list()

	for(j in 1:n.gene){
		Z = G[[j]]
		MAF = NULL
		for(i in 1:n.g){

			ID<-ID.position[[i]]
			Z1<-as.matrix(Z[ID,])
			MAF<- cbind(MAF, colMeans(Z1, na.rm = TRUE)/2 ) ##each column coresponds to a cohort, entries in the column are the MAFs of the variants in a gene region

		}
		MAF.list[[j]]<- MAF   ##each element in the MAF.list corresponds to a gene region
	}


	MAF.list.Groups<-list() ## each MAF.list.Groups[[g]] represent the ancestry-specific MAF for the g-th gene, 
					## each column in MAF.list.Groups[[g]] corresponds to a study
					## each row in MAF.list.Groups[[g]] corresponds to a variant the the g-th gene
	Map.Groups<-rep(0,n.g)
	

	# Get MAF.Groups when Group_Idx != NULL
	ID.Groups = unique(Group_Idx)
	n.ancestry = length(ID.Groups)	
	for(g in 1:n.gene){
		MAF.temp.Groups = NULL;
		for(j in 1:n.ancestry){

			temp<-which(Group_Idx == ID.Groups[j])
			Map.Groups[temp]<-j
			MAF.Groups = 0
			for(i in temp){
				MAF.Groups = MAF.Groups + MAF.list[[g]][,i] * n.each[i] / sum(n.each[temp])
			}
			MAF.temp.Groups = cbind(MAF.temp.Groups, MAF.Groups)
		}
		MAF.list.Groups[[g]] = MAF.temp.Groups
	}

	kernel = matrix(0, nrow = n.ancestry, ncol = n.ancestry)
	if(!is.Genetic.Similarity){
		diag(kernel) = 1
	}else{
		kernel = est.MAC(n.gene, n.ancestry, MAF.list.Groups)
	}

	return(kernel)


}
