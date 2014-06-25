#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.sph<-function(nmax,x){
	cn<-rep(0,nmax+1)  # Vector 
	rn<-rep(1,nmax+1)  # Vector
	cn[nmax+1]<-lcfe.sbl(nmax,x) # Last element
	rn[nmax+1]<-lcfe.sbd(nmax,x) # Last element
	Sn<-(0:(nmax+1))/x
	nj<-(nmax+1):2        # n+1
	rm<-rn
	cm<-cn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #rn[n-1]<-Sn[n+1]+cn[n]
   	  #cn[n-1]<-Sn[n-1]-1/rn[n-1]
   	  rn[n-1]<-Sn[n+1]+cn[n]
   	  cn[n-1]<-Sn[n-1]-1/rn[n-1]
   	  # modified (permits one step normalization)
   	  rm[n-1]<-Sn[n+1]*rm[n]+cm[n]
   	  cm[n-1]<-Sn[n-1]*rm[n-1]-rm[n]
   	  # Normalization
   	  #print(c(rn[n-1],cn[n-1]))
   	  if(abs(rm[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(rn[n-1],cn[n-1]))
   	  	 cm<-cm/rm[n-1] # this must be done first
   	  	 rm<-rm/rm[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(rn[n-1],cn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(rm[1])<abs(rm[2])){
      jn<-(rm/rm[1])*sin(x)/x # create functions for normalizations
   }else{
      jn<-(rm/rn[2])*(-cos(x)+sin(x)/x)*(1/x)   	
   }
   # Its Derivative
   if(abs(cm[1])>abs(cm[2])){
   	  djn<-(cm/cm[1])*(-jn[2])
   }else{
   	  djn<-(cm/cm[2])*(1/3)*(jn[1]-2*jn[3])
   }
   # Return results
   return(data.frame(rn,cn,jn,djn))
}
print(bess.sph(4,4))