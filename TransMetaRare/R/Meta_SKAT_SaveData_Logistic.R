Meta_SKAT_SaveData_Logistic <-
function(res, Z, X1, pi_1, res.out=NULL){

 
  # Get temp
  Q.Temp = t(res)%*%Z
  Q.Temp.Resampling<-NULL
  if(!is.null(res.out)){
 	Q.Temp.Resampling<-t(Z) %*% res.out
  }

  W.1 = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*%X1)%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1)) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  re<-list(Score=Q.Temp[1,], SMat.Summary = W.1,  MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)
}
