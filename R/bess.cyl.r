#-------------------------------------------------------------------------------
# CYLINDRICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
#---------------------------------------
# \gamma_{n-1}+\frac{1}{\gamma_n}=2S_{n}
# \gamma_{n-1}-\frac{gamma_n}=2D_{n}
#
# \gamma_{n}=S_{n+1}+D_{n+1}
# \frac{1}{\gamma_n}=S_n-D_n
#
# (S_n-D_n)(S_{n+1}+D_{n+1})=1
#---------------------------------------
bess.cyl<-function(nmax,x){ # PROBLEMAS COM ZEROS #
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(0,nmax+1)  # Vector
	Dn[nmax+1]<-lcfe.cbl(nmax,x)  # Last element - D_n=Dn[n+1]
	gn[nmax+1]<-lcfe.cbd(nmax,x)  # Last element - gamma_n=gn[n+1]
	Sn<-(0:nmax)/x                # S_n = Sn[n+1]
	nj<-(nmax+1):2                # n+1
	gm<-gn                        # gn modified - gm
	Dm<-Dn                        # Dn modified - Dm
	# DOWNWARD RECURRENCE         # Many choices
	RN<-1                         # Renormalization of modified functions counter
      for(n in nj){              # position (nmax+1):2 -> element nmax:1
   	   # gamma
         gn[n-1]=2*Sn[n]-1/gn[n]                 # [OK]
   	   Dn[n-1]<-Sn[n-1]-1/gn[n-1]              # [OK]
   	   # modified (permits one step normalization)
   	   gm[n-1]<-Sn[n]*gm[n]+Dm[n]
   	   Dm[n-1]<-Sn[n-1]*gm[n-1]-gm[n]
   	   # Normalization
   	   if(abs(gm[n-1])>1e100){
   	   	 cat("renorming...\n")
   	   	 #print(c(gn[n-1],Dn[n-1]))
   	   	 Dm<-Dm/gm[n-1] # this must be done first
   	   	 gm<-gm/gm[n-1] # otherwise the result will be wrong.
   	   }
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
   return(data.frame(gn,gm,Dn,Dm,besselJ(x,0:nmax)))
}
no<-4
xo<-5
cat("Dn    =",lcfe.cbl(no-1,xo),"\n")
cat("gamma =",lcfe.cbd(no-1,xo),"\n")
print(bess.cyl(no,xo))
