#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
sph.rjyn<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcf.sbld(nmax,x) # Last element
	gn[nmax+1]<-lcf.sbrd(nmax,x) # Last element
	Sn<-(0:(nmax+1))/x
	nj<-(nmax+1):2        # n+1
	Gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n+1]+Dn[n]
   	  #Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  Gm[n-1]<-Sn[n+1]+Dm[n]
   	  Dm[n-1]<-Sn[n-1]-1/Gm[n-1]
   	  # modified (permits one step normalization)
   	  gn[n-1]<-Sn[n+1]*gn[n]+Dn[n]
   	  Dn[n-1]<-Sn[n-1]*gn[n-1]-gn[n]
   	  # Normalization
   	  #print(c(gn[n-1],Dn[n-1]))
   	  if(abs(gn[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(gn[n-1],Dn[n-1]))
   	  	 Dn<-Dn/gn[n-1] # this must be done first
   	  	 gn<-gn/gn[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(gn[n-1],Dn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gn[1])<abs(gn[2])){
      jn<-(gn/gn[1])*sin(x)/x # create functions for normalizations
   }else{
      jn<-(gn/gn[2])*(-cos(x)+sin(x)/x)*(1/x)   	
   }
   # Its Derivative
   if(abs(Dn[1])>abs(Dn[2])){
   	  djn<-(Dn/Dn[1])*(-jn[2])
   }else{
   	  djn<-(Dn/Dn[2])*(1/3)*(jn[1]-2*jn[3])
   }
   # Return results
   return(data.frame(rh=Gm,An=Dm,jn,djn))
}
