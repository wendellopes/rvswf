#-------------------------------------------------------------------------------
# RICCATI BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.ric<-function(nmax,x,code="C"){
   if(abs(x)<1e-10){
      return(data.frame(Rn=rep(0,nmax+1),Dn=rep(0,nmax+1)))
   }
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("bess_ric",
            nmax=as.integer(nmax),
            x=as.double(x),
            Rn=as.double(dummy),
            Dn=as.double(dummy))
      return(data.frame(Rn=u$Rn,dRn=u$Dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   rm<-rn
	   Cn[nmax+1]<-lcfe.rbl(nmax,x) # Last element
	   rn[nmax+1]<-lcfe.sbd(nmax,x) # Last element
	   Cm<-Cn
	   # DOWNWARD RECURRENCE
	   RN<-1
      for(n in nmax:1){
      	  # original
      	  rn[n]<-S(n,x)+Cn[n+1]
      	  Cn[n]<-S(n,x)-1/rn[n]
      	  # modified (permits one step normalization)
      	  rm[n]<-S(n,x)*rm[n+1]+Cm[n+1]
      	  Cm[n]<-S(n,x)*rm[n]-rm[n+1]
      	  # Normalization
      	  if(abs(rm[n])>1e10){
      	  	 cat("renorming...\n")
      	  	 #print(c(rn[n-1],Cn[n-1]))
      	  	 Cm<-Cm/rm[n] # this must be done first
      	  	 rm<-rm/rm[n] # otherwise the result will be wrong.
      	  }
      	  #print(c(rn[n-1],Cn[n-1]))
      }
      # one step normalization taking care about zeros
      # Bessel function
      j0<-bess.zro(x)
      j1<-bess.uno(x)
      R0<-x*j0
      R1<-x*j1
      dR0<-cos(x)
      dR1<--j1+R0
      if(abs(R0)>1e-10){ # R_0 != 0
         jn<-(rm/rm[1])*R0 # create functions for normalizations
      }else{
         jn<-(rm/rm[2])*R1
      }
      # Its Derivative
      if(dR0>1e-10){ # R_0' != 0
      	  djn<-(Cm/Cm[1])*dR0
      }else{
      	  djn<-(Cm/Cm[2])*dR1
      }
      # Return results
      return(data.frame(Rn=jn,dRn=djn))
   }
}
