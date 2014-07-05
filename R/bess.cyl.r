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
   S<-function(n,x){
      return(n/x)
   }
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	gm<-gn                        # gn modified: gm -> J_{nmax}=1 -> J_{nmax}'=gn
	Dn[nmax+1]<-lcfe.cbl(nmax,x)  # Last element - D_n=Dn[n+1]
	gn[nmax+1]<-lcfe.cbd(nmax,x)  # Last element - gamma_n=gn[n+1]
	Dm<-Dn                        # Dn modified - Dm
	# DOWNWARD RECURRENCE         # Many choices
	RN<-1                         # Renormalization of modified functions counter
      for(n in nmax:1){              # position (nmax+1):2 -> element nmax:1
   	   # gamma
#         gn[n]=S(2*n,x)-1/gn[n+1]                 # [OK]
   	   gn[n]<-S(n,x)+Dn[n+1]
   	   Dn[n]<-S(n-1,x)-1/gn[n]              # [OK]
   	   # modified (permits one step normalization)
   	   gm[n]<-S(n,x)*gm[n+1]+Dm[n+1]
   	   Dm[n]<-S(n-1,x)*gm[n]-gm[n+1]
   	   # Normalization
   	   if(abs(gm[n])>1e10){
   	   	 cat("renorming...\n")
   	   	 Dm<-Dm/gm[n] # this must be done first
   	   	 gm<-gm/gm[n] # otherwise the result will be wrong.
   	   }
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gm[1])<abs(gm[2])){
      print("1")
      Jn<-(gm/gm[1])*besselJ(x,0) # create functions for normalizations
   }else{
      print("2")
      Jn<-(gm/gm[2])*besselJ(x,1)   	
   }
   # Its Derivative
   if(abs(Dm[1])>abs(Dm[2])){
      print("3")
   	  dJn<-(Dm/Dm[1])*(-Jn[2])
   }else{
      print("4")
   	  dJn<-(Dm/Dm[2])*.5*(Jn[1]-Jn[3])
   }
   # Return results
   return(data.frame(gn,Dn,Jn,dJn))
}
# no<-4
# xo<-5
# cat("Dn    =",lcfe.cbl(no-1,xo),"\n")
# cat("gamma =",lcfe.cbd(no-1,xo),"\n")
# print(bess.cyl(no,xo))
