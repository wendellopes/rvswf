#-------------------------------------------------------------------------------
# CYLINDRICAL WAVE GUIDE
#-------------------------------------------------------------------------------
vswf.cwg<-function(TM=TRUE,gama,kz,x,y,z,lmax,m,s=1){
   LMAX=lmax*(lmax+2)+1
   k<-sqrt(gama^2+kz^2)
   #----------------------------------------
   u<-vswf.qlm(kz/k,lmax)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   ll<-u$l
   mm<-u$m
   llp1<-1/sqrt(ll*(ll+1))
   llp1[1]<-0
   #----------------------------------------
   A<-2*(1i^ll)*((k/gama)^2)*Qlm*mm*llp1
   B<-2*(1i^(ll-1))*dQlm*llp1
   #----------------------------------------
   g<-2*pi*((-1i*s)^mm)*vswf.psi(gama,kz,x,y,z,m-s*mm,s)
   #----------------------------------------
   if(TM){# TM CWG
      GTE<- A*g
      GTM<--B*g
   }else{ # TE CWG
      GTE<-B*g
      GTM<-A*g
   }
   return(data.frame(GTE,GTM))
}
