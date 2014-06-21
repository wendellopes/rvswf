#-------------------------------------------------------------------------------
# Riccati Bessel Functions
#-------------------------------------------------------------------------------
tst.rbjn<-function(x,n){
	return(sqrt(pi/2)*besselJ(x,n+.5)*sqrt(x))
}