Meta_SKAT_SaveData_Linear <-
function(res, Z, X1, s2, res.out=NULL){

  # get Q
  Q.Temp = t(res)%*%Z

  Q.Temp.Resampling<-NULL
  if(!is.null(res.out)){
 	Q.Temp.Resampling<-t(Z) %*% res.out /s2
  }
  W.1 = t(Z) %*% Z - (t(Z) %*%X1)%*%solve(t(X1)%*%X1)%*% (t(X1) %*% Z ) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  
  #re<-list( Score=Q.Temp[1,] /sqrt(s2), SMat.Summary = W.1, MAF=MAF, Score.Resampling=Q.Temp.Resampling)
  re<-list( Score=Q.Temp[1,] /s2, SMat.Summary = W.1/s2, MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)

}
