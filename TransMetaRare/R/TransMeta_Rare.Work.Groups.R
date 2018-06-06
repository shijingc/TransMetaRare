TransMeta_Rare.Work.Groups <-
function(re, n.g, ID.Groups, Group_Idx){

	Score<-NULL
	SMat.Summary<-NULL
	Score.Resampling<-NULL
	
	n.Groups<-length(ID.Groups)	
	n.dim<-dim(re[[1]]$SMat.Summary)[1]
	n.dim.all<-n.Groups * n.dim
	SMat.Summary<-matrix(rep(0,n.dim.all*n.dim.all),ncol=n.dim.all)
	Score<-rep(0,n.dim.all)

	if(!is.null(re[[1]]$Score.Resampling)){
		temp1<-dim(re[[1]]$Score.Resampling)[2]
		Score.Resampling<-matrix(rep(0,n.dim.all*temp1),ncol=temp1)
	}

	idx.start=1
	for(i in 1:n.Groups){
		
		SMat.Summary1<-0
		Score1<-0
		Score1.Resampling<-0
		temp<-which(Group_Idx == ID.Groups[i])
		for(i in temp){
			SMat.Summary1 = SMat.Summary1 + re[[i]]$SMat.Summary
			Score1 = Score1 + re[[i]]$Score

			if(!is.null(Score.Resampling)){
				Score1.Resampling = Score1.Resampling + re[[i]]$Score.Resampling
				
			}
		}

		idx.end<-idx.start + n.dim -1
		idx<-idx.start:idx.end
		SMat.Summary[idx,idx]<-SMat.Summary1
		Score[idx]<-Score1	
		if(!is.null(Score.Resampling)){
			Score.Resampling[idx,]<-Score1.Resampling
		}
		idx.start = idx.end+1
	}

	re1<-list(Score=Score, SMat.Summary=SMat.Summary, Score.Resampling=Score.Resampling, m = n.dim, k = n.Groups)

	return(re1)

}
