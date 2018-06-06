matrix.power <-
function(A, n){
	a.eig = eigen(A, symmetric = TRUE)
	ID1 = which(a.eig$values > 0)
	a.power = a.eig$vectors[,ID1] %*% (a.eig$values[ID1]^n * t(a.eig$vectors[,ID1]))

	return(a.power)
}
