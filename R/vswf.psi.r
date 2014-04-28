#-------------------------------------------------------------------------------
# PSI(k,r)
#-------------------------------------------------------------------------------
vswf.psi<-function(gama,kz,x,y,z,m,s=1){
   rho<-sqrt(x^2+y^2)
   if((x==0)&&(y==0)){
      cph<-1
      sph<-0
   }else{
      cph<-x/rho
      sph<-y/rho   
   }
   u<-besselJ(gama*rho,m)*((cph+1i*s*sph)^m)*exp(1i*kz*z)
   return(u)
}
