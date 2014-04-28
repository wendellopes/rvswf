#-------------------------------------------------------------------------------
# GENERIC PLANE WAVE
#-------------------------------------------------------------------------------
vswf.gpw<-function(kx,ky,kz,ux,uy,uz,lmax,xo=0,yo=0,zo=0){
   a<-vswf.vsh(kx,ky,kz,lmax) # Calculo direto de X_{lm}
   k<-sqrt(Conj(kx)*kx+Conj(ky)*ky+Conj(kz)*kz) # k real
   u<-sqrt(Conj(ux)*ux+Conj(uy)*uy+Conj(uz)*uz) # u imag
   LMAX=lmax*(lmax+2)+1
   # k normalization
   hkx<-kx/k
   hky<-ky/k
   hkp<-(hkx+1i*hky)/sqrt(2)
   hkm<-(hkx-1i*hky)/sqrt(2)
   hkz<-kz/k
   # u normalization
   hux<-ux/u
   huy<-uy/u
   huz<-uz/u
   hup<-(hux+1i*huy)/sqrt(2)
   hum<-(hux-1i*huy)/sqrt(2)
   # ALFA a = k x u
   ham<-1i*(hkm*huz-hkz*hum)
   haz<-1i*(hkp*hum-hkm*hup)
   hap<-1i*(hkz*hup-hkp*huz)
   # Normalized Associated Legendre Polynomials
   U<-vswf.qlm(kz/k,lmax+1)
   qlm<-U$Qlm
   rm(U)
   # Constants
   cplm<-function(l,m){return(sqrt(l*(l+1)-m*(m+1)))}
   cmlm<-function(l,m){return(sqrt(l*(l+1)-m*(m-1)))}
   il<-function(l){return((1i^l)/sqrt(l*(l+1)))}
   # Calculations
   GTE.ex<-rep(0,LMAX)
   GTM.ex<-rep(0,LMAX)
   Il<-rep(0,LMAX)
   for(l in 1:lmax){
      m<--l:l
      lmm<-vswf.jlm(l,m-1)
      lmz<-vswf.jlm(l,m)
      lmp<-vswf.jlm(l,m+1)
      Il[lmz]<-1i^l
      if((kx==0)&&(ky==0)){
         em<-1
         ez<-1
         ep<-1
      }else{
         gk<-sqrt(kx^2+ky^2)
         czt<-kx/gk
         szt<-ky/gk
         em<-(czt-1i*szt)^(m-1)
         ez<-(czt-1i*szt)^m
         ep<-(czt-1i*szt)^(m+1)
      }
      # Explicit calculation
      GTE.ex[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hup/sqrt(2)+
         m*ez*qlm[lmz]*huz+cmlm(l,m)*em*qlm[lmm]*hum/sqrt(2))
      GTM.ex[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hap/sqrt(2)+
         m*ez*qlm[lmz]*haz+cmlm(l,m)*em*qlm[lmm]*ham/sqrt(2))
   } 
   # Implicit calculation
   GTE.in<-4*pi*Il*(Conj(a$Xm)*hum+Conj(a$Xz)*huz+Conj(a$Xp)*hup)
   GTM.in<-4*pi*Il*(Conj(a$Xm)*ham+Conj(a$Xz)*haz+Conj(a$Xp)*hap)
   return(data.frame(GTE.ex,GTM.ex,GTE.in,GTM.in))
}
