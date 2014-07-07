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
lcfa.cyl<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("lcfa_cyl",
            nmax=as.integer(nmax),
            x=as.double(x),
            gn=as.double(dummy),
            Dn=as.double(dummy),
            NMAX=as.integer(NMAX))
      return(data.frame(gn=u$gn,Dn=u$Dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Dn<-rep(0,nmax+1)  # Vector 
	   gn<-rep(1,nmax+1)  # Vector
	   Dn[nmax+1]<-lcfe.cbl(nmax,x)  # Last element - D_n=Dn[n+1]
	   gn[nmax+1]<-lcfe.cbd(nmax,x)  # Last element - gamma_n=gn[n+1]
	   # DOWNWARD RECURRENCE         # Many choices
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
#        gn[n]=S(2*n,x)-1/gn[n+1]   # [OK]
      	gn[n]<-S(n,x)+Dn[n+1]      # [OK]
      	Dn[n]<-S(n-1,x)-1/gn[n]    # [OK]
      }
      return(data.frame(gn,Dn))
   }
}
# no<-4
# xo<-5
# cat("Dn    =",lcfe.cbl(no-1,xo),"\n")
# cat("gamma =",lcfe.cbd(no-1,xo),"\n")
# print(bess.cyl(no,xo))
