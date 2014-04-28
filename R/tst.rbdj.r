#-------------------------------------------------------------------------------
# Derivative of Riccati Bessel Functions
#-------------------------------------------------------------------------------
tst.rbdj<-function(x,n){
   a<-n/(2*n+1)
   b<-(n+1)/(2*n+1)
	return(b*tst.sbjn(x,n-1)-a*tst.sbjn(x,n+1))
}