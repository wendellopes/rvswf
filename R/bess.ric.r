#-------------------------------------------------------------------------------
# RICCATI BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.ric<-function(nmax,x){
	Cn<-rep(0,nmax+1)  # Vector 
	rn<-rep(1,nmax+1)  # Vector
	Cn[nmax+1]<-lcfe.rbl(nmax,x) # Last element
	rn[nmax+1]<-lcfe.sbd(nmax,x) # Last element
	Sn<-(0:(nmax))/x
	nj<-(nmax+1):2        # n+1
	rm<-rn
	Cm<-Cn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #rn[n-1]<-Sn[n]+Cn[n]
   	  #Cn[n-1]<-Sn[n]-1/rn[n-1]
   	  rn[n-1]<-Sn[n]+Cn[n]
   	  Cn[n-1]<-Sn[n]-1/rn[n-1]
   	  # modified (permits one step normalization)
   	  rm[n-1]<-Sn[n]*rm[n]+Cm[n]
   	  Cm[n-1]<-Sn[n]*rm[n-1]-rm[n]
   	  # Normalization
   	  #print(c(rn[n-1],Cn[n-1]))
   	  if(abs(rm[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(rn[n-1],Cn[n-1]))
   	  	 Cm<-Cm/rm[n-1] # this must be done first
   	  	 rm<-rm/rm[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(rn[n-1],Cn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(rm[1])<abs(rm[2])){
      jn<-(rm/rm[1])*sin(x) # create functions for normalizations
   }else{
      jn<-(rm/rm[2])*(sin(x)/x-cos(x)) 	
   }
   # Its Derivative
   if(abs(Cm[1])>abs(Cm[2])){
   	  djn<-(Cm/Cm[1])*(-jn[2])
   }else{
   	  djn<-(Cm/Cm[2])*(1/3)*(2*jn[1]-jn[3])
   }
   # Return results
   return(data.frame(rn,Cn,jn,djn))
}
