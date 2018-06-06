Meta_SKAT_SaveData_Kinship <-
function(res, Z, P1){


  # get Q
  Q.Temp = t(res)%*%Z

  Q.Temp.Resampling<-NULL
  W.1 = t(Z) %*% (P1 %*% Z) # t(Z) P0 Z
  MAF = colMeans(Z)/2

  re<-list( Score=Q.Temp[1,] , SMat.Summary = W.1, MAF=MAF, Score.Resampling=Q.Temp.Resampling)  
  return(re)


}
