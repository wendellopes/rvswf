#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
#---------------------------------------
# \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}
# n\rho_{n-1}-(n+1)\frac{rho_n}=(2n+1)c_{n}
#
# \rho_{n}=S_{n+2}+c_{n+1}
# \frac{1}{\rho_n}=S_n-c_n
#
# (S_n-c_n)(S_{n+2}+c_{n+1})=1
#---------------------------------------
lcfa.sph<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("lcfa_sph",
            nmax=as.integer(nmax),
            x=as.double(x),
            rn=as.double(dummy),
            cn=as.double(dummy),
            NMAX=as.integer(NMAX))
      return(data.frame(rn=u$rn,cn=u$cn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   cn[nmax+1]<-lcfe.sbl(nmax,x)  # Last element - c_n=cn[n+1]
	   rn[nmax+1]<-lcfe.sbd(nmax,x)  # Last element - rho_n=rn[n+1]
	   # DOWNWARD RECURRENCE     
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
      	rn[n]<-S(n+1,x)+cn[n+1]    # [OK]
      	cn[n]<-S(n-1,x)-1/rn[n]    # [OK]
      }
      return(data.frame(rn,cn))
   }
}
