#-------------------------------------------------------------------------------
# CYLINDRICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.cyl<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcf.cbld(nmax,x) # Last element
	gn[nmax+1]<-lcf.cbrd(nmax,x) # Last element
	Sn<-(0:nmax)/x
	nj<-(nmax+1):2        # n+1
	Gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n]+Dn[n]
   	  #Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  Gm[n-1]<-Sn[n]+Dm[n]
   	  Dm[n-1]<-Sn[n-1]-1/Gm[n-1]
   	  # modified (permits one step normalization)
   	  gn[n-1]<-Sn[n]*gn[n]+Dn[n]
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
      Jn<-(gn/gn[1])*besselJ(x,0) # create functions for normalizations
   }else{
      Jn<-(gn/gn[2])*besselJ(x,1)   	
   }
   # Its Derivative
   if(abs(Dn[1])>abs(Dn[2])){
   	  dJn<-(Dn/Dn[1])*(-Jn[2])
   }else{
   	  dJn<-(Dn/Dn[2])*.5*(Jn[1]-Jn[3])
   }
   #DIRECT CALCULATION
   Jl<-dJl<-Jn
   Jl[1]<-besselJ(x,0)
   dJl[1]<-Jl[1]*Dm[1]
   for(n in 1:(nmax)){
   	  Jl[n+1]<-Jl[n]/Gm[n]
   	  dJl[n+1]<-Jl[n+1]*Dm[n+1]
   }
   # Return results
   return(data.frame(gm=Gm,Dn=Dm,Jn,Jl,dJn,dJl))
}
