#-------------------------------------------------------------------------------
# BESSEL BEAM Z
#-------------------------------------------------------------------------------
vswf.bbz<-function(gama,kz,xo,yo,zo,lmax,M,s){
   k<-sqrt(gama^2+kz^2)
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-vswf.qlm(kz/k,lmax+1)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   rm(u)
   #----------------------------------------
   GTE<-rep(0,LMAX)
   GTM<-rep(0,LMAX)
   for(l in 1:lmax){
      m<--l:l
      psio<-vswf.psi(gama,kz,xo,yo,zo,M-s*m,s)
      if(l>0){
         A0<-4*pi*(1i^(l-s*m))*psio/sqrt(l*(l+1))
      }else{
         A0<-0
      }
     GTE[vswf.jlm(l,m)]<-m*Qlm[vswf.jlm(l,m)]*A0
     GTM[vswf.jlm(l,m)]<-1i*((gama/k)^2)*dQlm[vswf.jlm(l,m)]*A0
   }
   return(data.frame(GTE,GTM))
}
