#-------------------------------------------------------------------------------
# CYLINDRICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.cyl<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcfe.cbl(nmax,x) # Last element
	gn[nmax+1]<-lcfe.cbd(nmax,x) # Last element
	Sn<-(0:nmax)/x
	nj<-(nmax+1):2        # n+1
	gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n]+Dn[n]
   	  #Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  gn[n-1]<-Sn[n]+Dn[n]
   	  Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  # modified (permits one step normalization)
   	  gm[n-1]<-Sn[n]*gm[n]+Dm[n]
   	  Dm[n-1]<-Sn[n-1]*gm[n-1]-gm[n]
   	  # Normalization
   	  #print(c(gn[n-1],Dn[n-1]))
   	  if(abs(gm[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(gn[n-1],Dn[n-1]))
   	  	 Dm<-Dm/gm[n-1] # this must be done first
   	  	 gm<-gm/gm[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(gn[n-1],Dn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gm[1])<abs(gm[2])){
      Jn<-(gm/gm[1])*besselJ(x,0) # create functions for normalizations
   }else{
      Jn<-(gm/gm[2])*besselJ(x,1)   	
   }
   # Its Derivative
   if(abs(Dm[1])>abs(Dm[2])){
   	  dJn<-(Dm/Dm[1])*(-Jn[2])
   }else{
   	  dJn<-(Dm/Dm[2])*.5*(Jn[1]-Jn[3])
   }
   # Return results
   return(data.frame(gn,Dn,Jn,dJn))
}
