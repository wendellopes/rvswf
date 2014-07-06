#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
bess.zro<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(x>1){
      return(sin(x)/x)
   }else{
      S<-1
      for(j in seq(10,1,-1)){
         S<-1-Kj(x,2*j)*Kj(x,2*j+1)*S
      }
      return(S)
   }
}
#-------------------------------------------------------------------------------
bess.uno<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(x>1){
      return((sin(x)/x-cos(x))/x)
   }else{
      S<-0
      for(j in seq(10,1,-1)){
         S<-(2*j)/(2*j+1)-Kj(x,2*j+1)*Kj(x,2*j+2)*S
      }
      S<-x*S/2
      return(S)
   }
}
#-------------------------------------------------------------------------------
bess.sph<-function(nmax,x,code="C"){
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("bess_sph",
            nmax=as.integer(nmax),
            x=as.double(x),
            jn=as.double(dummy),
            dn=as.double(dummy))
      return(data.frame(jn=u$jn,djn=u$dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   rm<-rn
	   cn[nmax+1]<-lcfe.sbl(nmax,x) # Last element
	   rn[nmax+1]<-lcfe.sbd(nmax,x) # Last element
	   cm<-cn
	   # DOWNWARD RECURRENCE
	   RN<-1
      for(n in nmax:1){
         #print(c(n,rm[n+1],cm[n+1]))
         #rn[n]<-S(2*n+1,x)-1/rn[n+1]  # Only \rho_n
         rn[n]<-S(n+1,x)+cn[n+1]  # [OK]
         cn[n]<-S(n-1,x)-1/rn[n]  # [OK]
      	# modified (permits one step normalization)
         # AQUI! #
      	rm[n]<-S(n+1,x)*rm[n+1]+cm[n+1]
      	cm[n]<-S(n-1,x)*rm[n]-rm[n+1]
         # zeroing
         if(abs(rm[n])<1e-14){
            rm[n]<-0
         }
         if(abs(cm[n])<1e-14){
            cm[n]<-0
         }
         #print(c(n-1,rm[n],cm[n]))
      	# Normalization
      	#print(c(rn[n-1],cn[n-1]))
      	if(abs(rm[n])>1e2){
      	  	#cat("renorming...\n")
      	  	#print(c(rn[n-1],cn[n-1]))
      	  	cm<-cm/rm[n] # this must be done first
      	   rm<-rm/rm[n] # otherwise the result will be wrong.
      	}
      }
      # one step normalization taking care about zeros
      # Bessel function
      b0<-bess.zro(x)
      b1<-bess.uno(x)
      if(abs(b0)>1e-10){ # j_0 != 0
         jn<-(rm/rm[1])*b0 # create functions for normalizations
      }else{            # j_0 == 0
         jn<-(rm/rm[2])*b1
      }
      # Its Derivative
      if(abs(b1)>1e-10){ # j_0' != 0
      	  djn<-(cm/cm[1])*(-b1) # j_0'=-j_1
      }else{            # j_0' == 0
      	  djn<-(cm/cm[2])*b0 # j_1' = (2/x)j_0'+j_0 = j_0 # APROXIMADO!
      }
      return(data.frame(rn,cn,jn,djn))
   }
   # Return results
}
