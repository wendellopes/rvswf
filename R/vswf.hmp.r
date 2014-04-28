#-------------------------------------------------------------------------------
# MULIPOLOS DE HANSEN
#-------------------------------------------------------------------------------
vswf.hmp<-function(k,x,y,z,lmax){
#------------------------------------------------
   LMAX=lmax*(lmax+2)+1
   u<-vswf.vsh(x,y,z,lmax) 
   r<-sqrt(x^2+y^2+z^2)
#------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS CALCULATIONS
#------------------------------------------------
# GNU SCIENTIFIC LIBRARY (GSL) - RESULT AS MATRIX
#   library(gsl)
#   jl<-bessel_jl_steed_array(lmax+1,k*r) # Funcoes de bessel de 0 a lmax+1
#------------------------------------------------
#LENTZ CONTINUED FRACTION AND DOWNWARD RECURRENCE
#source("vswf.sbf.R")
jll<-vswf.sbf(lmax+1,k*r)
#jl<-jll$jn # FUNCIONA TB, RESULT AS NUMERIC 
jl<-as.matrix(jll$jn) # RESULT AS MATRIX
#------------------------------------------------
   jl.m<-0 # Correct Value: cos(k*r)/(k*r); Nevertheless, M,N,L starts in 1
   jl.o<-jl[1]
   jl.M<-jl[2]
#------------------------------------------------
   for(l in 1:lmax){
      n.rep<-2*l+1
      jl.m<-c(jl.m,rep(jl[l  ],n.rep))
      jl.o<-c(jl.o,rep(jl[l+1],n.rep))
      jl.M<-c(jl.M,rep(jl[l+2],n.rep))
   }
#------------------------------------------------
   # Constants
   K<-2*u$l+1
   Kl<-u$l/K
   KL<-(u$l+1)/K
   Kll<-sqrt(u$l*(u$l+1))
   kl<-sqrt(Kl)
   kL<-sqrt(KL)
   p1<-KL*jl.m-Kl*jl.M # Derivative of jl(x)
   p2<-(jl.m+jl.M)/K   # jl(x)/x
#------------------------------------------------
   # M
   M.m<-jl.o*u$Xm
   M.z<-jl.o*u$Xz
   M.p<-jl.o*u$Xp
#------------------------------------------------
   # N
   N.m<-p1*u$Vm+Kll*p2*u$Ym
   N.z<-p1*u$Vz+Kll*p2*u$Yz
   N.p<-p1*u$Vp+Kll*p2*u$Yp
#------------------------------------------------
#   # N_{lm} by Y_{j,l}^m (shorter expression)
#   Ny.m<-kL*jl.m*u$Ym.lm-kl*jl.M*u$Ym.Lm
#   Ny.z<-kL*jl.m*u$Yz.lm-kl*jl.M*u$Yz.Lm
#   Ny.p<-kL*jl.m*u$Yp.lm-kl*jl.M*u$Yp.Lm
#------------------------------------------------
   U<-data.frame(M.m,M.z,M.p,N.m,N.z,N.p)
#                 ,Ny.m,Ny.z,Ny.p)
}
