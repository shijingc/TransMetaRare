SKAT_2Kernel_Copula_Optimal_PValue <-
function(u, Sigma){


	obj.norm.Copula = normalCopula(Sigma[lower.tri(Sigma)],dim = nrow(Sigma),dispstr = "un")
  	#obj.t.Copula = tCopula(Sigma[lower.tri(Sigma)],dim = nrow(Sigma),dispstr = "un")
    	minP <- min(u)
  
  	pvalue.normCopula <- 1- pCopula(1-minP, obj.norm.Copula)
	#pvalue.tCopula <- 1- pCopula(1-minP, obj.t.Copula)
	
	#res = c(pvalue.normCopula, pvalue.tCopula )
	res = pvalue.normCopula
	return(res)

}
