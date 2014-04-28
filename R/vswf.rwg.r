#-------------------------------------------------------------------------------
# RECTANGULAR WAVE GUIDE
#-------------------------------------------------------------------------------
vswf.rwg<-function(TM=TRUE,kx,ky,kz,x,y,z,lmax){
   if(TM){
      s<--1   #s=-1 -> TM MODE
   }else{
      s<- 1   #s=+1 -> TE MODE
   }
   LMAX=lmax*(lmax+2)+1
   gama<-sqrt(kx^2+ky^2)
   k<-sqrt(kx^2+ky^2+kz^2)
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
   EXmY<-exp(1i*(kx*x-ky*y))
   EXpY<-exp(1i*(kx*x+ky*y))
   czt<-kx/gama
   szt<-ky/gama
   eimz<-(czt+1i*szt)^mm    
   f<-(-1)^mm
   #----------------------------------------
   g<-pi*exp(1i*kz*z)*((EXmY+f*Conj(EXmY))*eimz+
      s*((EXpY+f*Conj(EXpY))*Conj(eimz)))/2
   if(TM){# TM CWG
      GTE<- A*g
      GTM<--B*g
   }else{ # TE CWG
      GTE<-B*g
      GTM<-A*g
   }
   return(data.frame(GTE,GTM))
}
