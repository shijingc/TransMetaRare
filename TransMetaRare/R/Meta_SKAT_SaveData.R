Meta_SKAT_SaveData <-
function(Z, obj.res, SetID, impute.method = "fixed"){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	re<-1

	out.z<-Meta_SKAT_MAIN_Check_Z(Z, n, obj.res$id_include, SetID, impute.method=impute.method)
	
	if(class(obj.res)== "SKAT_NULL_Model_EMMAX"){
		out = Meta_SKAT_SaveData_Kinship(obj.res$res, out.z$Z.test, obj.res$P)

	} else if(obj.res$out_type == "C"){
		out = Meta_SKAT_SaveData_Linear(obj.res$res,out.z$Z.test
			,obj.res$X1, obj.res$s2, obj.res$res.out)
		  
	} else if (obj.res$out_type == "D"){

		out = Meta_SKAT_SaveData_Logistic(obj.res$res, out.z$Z.test
			,obj.res$X1, obj.res$pi_1, obj.res$res.out)
		
	}

	re=list(Score=out$Score, SMat.Summary = out$SMat.Summary, MAF=out.z$MAF, missing_rate=out.z$missing_rate,
	Score.Resampling=out$Score.Resampling )
	
	return(re)

}
