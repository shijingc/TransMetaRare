TransMeta_Rare_Null_Model_EmmaX <-
function(y.list, x.list, n.cohort, out_type="C", K.list=NULL, Kin.File.list=NULL){


	re<-list()
	re$n.g=n.cohort
	re$out<-list()
	re$ID<-list()
	re$n.each<-rep(0,n.cohort)
	re$out_type=out_type
	count<-1
	n.x = length(x.list)
	for(i in 1:n.cohort){
		y = y.list[[i]]
		X = x.list[[i]]
		if(is.vector(X)){
			if(X[1] == "intercept"){
				formula = y ~ 1
			} else if( length(X) == 1 && is.character(X)){
				msg<-sprintf("ERROR: %d element of x.list is wrong. For intercept-only, it should be \"intercept\".",i)
				stop(msg)
			} else {
				formula = y ~ X
			}
		} else {
			formula = y ~ X
		}
		
		if(!is.null(K.list)){
				re$out[[i]]<-SKAT::SKAT_NULL_emmaX(formula, K=K.list[[i]])
			} else if(!is.null(Kin.File.list)){
				re$out[[i]]<-SKAT::SKAT_NULL_emmaX(formula, Kin.File=Kin.File.list[[i]])
			} else {
				stop("Either K.list or Kin.File.list should have non-null value!")
			}


		if(length(re$out[[i]]$id_include) != length(y)){
			msg<-sprintf("Missing values in cohort %d. The current version cannot handle missing values in either phenotypes or covariates!",i)
			stop(msg)
		}

		n1<-length(re$out[[i]]$res)
		re$ID[[i]]<-count:(count +n1-1)
		count<-count + n1

		re$n.each[i]<-n1
	}

	class(re)<-"META_NULL_Model_EmmaX"



	if(out_type == "D"){
		dd = list()
		x.names = NULL
		for(j in 1:dim(x.list[[1]])[2]){
			x.names = c(x.names, paste("x", j, sep ='') )
		}

		for(i in 1:n.g){
			dd[[i]] = data.frame(y.list[[i]], x.list[[i]])
			colnames(dd[[i]]) = c('y', x.names) 
			if(!is.null(K.list)){

				H0 = glmmkin(y ~ ., data = dd[[i]],  kins = K.list[[i]], family = binomial(link = 'logit'))
				re$out[[i]]$res <- H0$residuals
				re$out[[i]]$P <- H0$P
   
			} else if(!is.null(Kin.File.list)){

				H0 = glmmkin(y ~ ., data = dd[[i]],  kins = K.list[[i]], family = binomial(link = 'logit'))
				re$out[[i]]$res <- H0$residuals
				re$out[[i]]$P <- H0$P

			}
		}

	}


	return(re)
	
}
