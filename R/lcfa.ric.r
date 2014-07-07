#-------------------------------------------------------------------------------
# RICCATI-BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
#---------------------------------------
# \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}
# (n+1)\rho_{n-1}-n\frac{1}{rho_n}=(2n+1)C_{n}
#
# \rho_{n}=S_{n+1}+C_{n+1}
# \frac{1}{\rho_n}=S_{n+1}-C_n
#
# (S_{n+1}-C_n)(S_{n+1}+C_{n+1})=1
#---------------------------------------
lcfa.ric<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("lcfa_ric",
            nmax=as.integer(nmax),
            x=as.double(x),
            rn=as.double(dummy),
            Cn=as.double(dummy),
            NMAX=as.integer(NMAX))
      return(data.frame(rn=u$rn,Dn=u$Cn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   Cn[nmax+1]<-lcfe.rbl(nmax,x)  # Last element - c_n=Cn[n+1]
	   rn[nmax+1]<-lcfe.sbd(nmax,x)  # Last element - rho_n=rn[n+1]
	   # DOWNWARD RECURRENCE     
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
      	rn[n]<-S(n,x)+Cn[n+1]      # [OK]
      	Cn[n]<-S(n,x)-1/rn[n]      # [OK]
      }
      return(data.frame(rn,Dn=Cn))
   }
}
