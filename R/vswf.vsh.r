#' Calculates the Vector Spherical Harmonics from zero to lmax.
#' 
#' @details Calculate the three complex components \eqn{\hat{e}_-},
#' \eqn{\hat{e}_z} and \eqn{\hat{e}_+} of the three vector spherical
#' harmonics \eqn{\bm{X}_{lm}}, \eqn{\bm{V}_{lm}} and \eqn{\bm{Y}_{lm}} and
#' also the three complex components of the three vector spherical harmonics
#' \eqn{\bm{Y}_{l,l}^m}, \eqn{\bm{Y}_{l,l-1}^m} and \eqn{\bm{Y}_{l,l+1}^m}.
#' @param x The component \eqn{x} of the position vector.
#' @param y The component \eqn{y} of the position vector.
#' @param z The component \eqn{z} of the position vector.
#' @param lmax The maximum value of \eqn{l}.
#' @return A list with the values of \eqn{l}, \eqn{-l\leq m\geq l}, the scalar
#' spherical harmonics and the components of the vector spherical harmonics.
#' @include vswf.jlm.r
#' @export
#' @examples
#' th<-pi/3
#' ph<-pi/4
#' VSH<-vswf.vsh(sin(th)*cos(ph),sin(th)*sin(ph),cos(th),5)
#' str(VSH)
vswf.vsh<-function(x,y,z,lmax,code="C"){
   LMAX=lmax*(lmax+2)+1
   dummy<-rep(0,LMAX)
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      if(TM){
         tm<-1
      }else{
         tm<-0
     }
      u<-.C("vswf_vsh",
         x=as.double(x),
         y=as.double(y),
         z=as.double(z),

         lmax=as.integer(lmax),

         Y=as.complex(dummy),

         VM=as.complex(dummy),
         VZ=as.complex(dummy),
         VP=as.complex(dummy),

         XM=as.complex(dummy),
         XZ=as.complex(dummy),
         XP=as.complex(dummy),

         YM=as.complex(dummy),
         YZ=as.complex(dummy),
         YP=as.complex(dummy))
      return(data.frame(Y=u$Y,VM=u$VM,VZ=u$VZ,VP=u$VP,
                              XM=u$XM,XZ=u$XZ,XP=u$XP,
                              YM=u$YM,YZ=u$YZ,YP=u$YP))
   }else{
#-------------------------------------------------------------------------------
# UNIDADES BASICAS
#-------------------------------------------------------------------------------
      rho<-sqrt(x^2+y^2)  # rho
      r<-sqrt(rho^2+z^2)  # r
      cth<-z/r            # cos(theta)
      sth<-rho/r          # sin(theta)
      cph<-x/rho          # cos(phi)
      sph<-y/rho          # sin(phi)
      if((x==0)&&(y==0)){
         cph<-1
         sph<-0
      }
#-------------------------------------------------------------------------------
# Constants for Qlm
#-------------------------------------------------------------------------------
      alfaQ<-function(l,m){
         return(sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m))))
      }
      betaQ<-function(l,m){
         return(sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m))))
      }
      gammaQ<-function(l){
         return(sqrt((2*l+1)/(2*l)))
      }
      deltaQ<-function(l){
         return(sqrt(2*l+1))
      }
#-------------------------------------------------------------------------------
# Constants for X
#-------------------------------------------------------------------------------
      cp<-function(l,m){
         return(sqrt(l*(l+1)-m*(m+1)))
      }
      cm<-function(l,m){
         return(sqrt(l*(l+1)-m*(m-1)))
      }
      cz<-function(l,m){
         return(m)
      }
      cl<-function(l,m){
         return(l)
      }
      ll<-function(l,m){
         return(1/sqrt(l*(l+1)))
      }
#-------------------------------------------------------------------------------
# Constants for Y e V
#-------------------------------------------------------------------------------
      Km<-function(l,m,q){
         return(sqrt((l-m)*(l-q))/sqrt((2*l-1)*(2*l+1)))
      }
      Ko<-function(l,m,q){
         return(sqrt((l-m)*(l+q))/sqrt((2*l-1)*(2*l+1)))
      }
      Kp<-function(l,m,q){
         return(sqrt((l+m)*(l+q))/sqrt((2*l-1)*(2*l+1)))
      }
#-------------------------------------------------------------------------------
# Constants for Y.lm e Y.LM (Blatt and Weisskopf)
#-------------------------------------------------------------------------------
      KYmlm<-function(l,m){
         return(sqrt((l+m-1)*(l+m))/sqrt(2*l*(2*l-1)))
      }
      KYolm<-function(l,m){
         return(sqrt((l-m)*(l+m))/sqrt(l*(2*l-1)))
      }
      KYplm<-function(l,m){
         return(sqrt((l-m-1)*(l-m))/sqrt(2*l*(2*l-1)))
      }
      #
      KYmLm<-function(l,m){
         return(sqrt((l-m+1)*(l-m+2))/sqrt(2*(l+1)*(2*l+3)))
      }
      KYoLm<-function(l,m){
         return(sqrt((l-m+1)*(l+m+1))/sqrt((l+1)*(2*l+3)))
      }
      KYpLm<-function(l,m){
         return(sqrt((l+m+2)*(l+m+1))/sqrt(2*(l+1)*(2*l+3)))
      }
#-------------------------------------------------------------------------------
# First scalar spherical harmonics
#-------------------------------------------------------------------------------
      if(lmax<1){lmax<-1}                         # At least 1 term
      LMAX=lmax*(lmax+2)+1                        # Vector for lmax
      LMAXE=(lmax+1)*(lmax+3)+1                   # Vector (lmax+1) - Extended
#-------------------------------------------------------------------------------
      Qlm<-rep(0,LMAXE)                           # First 4 terms for Qlm
      Qlm[vswf.jlm(0,0)]<-1/sqrt(4*pi)                      # Q00
      Qlm[vswf.jlm(1,1)]<--gammaQ(1)*sth*Qlm[vswf.jlm(0,0)] # Q11
      Qlm[vswf.jlm(1,0)]<-sqrt(3)*cth*Qlm[1]                # Q10
      Qlm[vswf.jlm(1,-1)]<--Qlm[vswf.jlm(1,1)]              # Q11*(-1)
#-------------------------------------------------------------------------------
      Ylm<-rep(1,LMAXE)
      eimphi.m<-rep(1,lmax+1)                 # m=0
      eimphi.p<-rep(1,lmax+1)
      eimphi.m[2]<-eimphi.m[1]*(cph-1i*sph)   # m=1
      eimphi.p[2]<-eimphi.p[1]*(cph+1i*sph)
#-------------------------------------------------------------------------------
      Ylm[1:4]<-Qlm[1:4]*c(1,(cph-1i*sph),1,(cph+1i*sph))
#-------------------------------------------------------------------------------
# VECTORS 
#-------------------------------------------------------------------------------
      Y.l.m<-dummy
      Y.l.o<-dummy
      Y.l.M<-dummy
#------------------
      Y.o.m<-dummy
      Y.o.o<-dummy
      Y.o.M<-dummy
#------------------
      Y.L.m<-dummy
      Y.L.o<-dummy
      Y.L.M<-dummy
#------------------
# VECTOR SPH. HARM.
#------------------
      Vm<-dummy
      Vz<-dummy
      Vp<-dummy
#------------------
      Xm<-dummy
      Xz<-dummy
      Xp<-dummy
#------------------
      Ym<-dummy
      Yz<-dummy
      Yp<-dummy
#------------------
# CONSTANTES
#------------------
      c.l<-dummy
      l.l<-dummy
#------------------
      c.p<-dummy
      c.m<-dummy
      c.z<-dummy
#------------------
      KmLm<-dummy
      KoLo<-dummy
      KpLM<-dummy
#------------------
      KmlM<-dummy
      Kolo<-dummy
      Kplm<-dummy
#------------------
      Km.lm<-dummy
      Ko.lm<-dummy
      Kp.lm<-dummy
#------------------
      Km.Lm<-dummy
      Ko.Lm<-dummy
      Kp.Lm<-dummy
#-------------------------------------------------------------------------------
# KERNEL --- CALCULATION OF VECTORS
#-------------------------------------------------------------------------------
      for(l in 1:lmax){
         mm<-(-l):l 
#-------------------------------------------------------------------------------
         # NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS & SCALAR SPHERICAL HARMONICS
         l<-l+1                              # l stars in 2 (deg 5)
         m.p<-0:(l-2)
         m.m<-(-l):(-1)
         csp<-(-1)^m.m
         Qlm[vswf.jlm(l,l  )]=-gammaQ(l)*sth*Qlm[vswf.jlm(l-1,l-1)]     #OK
         Qlm[vswf.jlm(l,l-1)]=deltaQ(l)*cth*Qlm[vswf.jlm(l-1,l-1)]      #OK
         Qlm[vswf.jlm(l,m.p)]=alfaQ(l,m.p)*cth*Qlm[vswf.jlm(l-1,m.p)]-
            betaQ(l,m.p)*Qlm[vswf.jlm(l-2,m.p)]                         #OK
         Qlm[vswf.jlm(l,m.m)]=csp*Qlm[vswf.jlm(l,abs(m.m))]             #OK
         # EXP(I*L*PHI)
         eimphi.p[l+1]<-eimphi.p[l]*(cph+1i*sph)
         eimphi.m[l+1]<-eimphi.m[l]*(cph-1i*sph)
         eimphi<-c(rev(eimphi.m[2:(l+1)]),eimphi.p[1:(l+1)])
         # Ylm
         Ylm[vswf.jlm(l,-l):vswf.jlm(l,l)]<-Qlm[vswf.jlm(l,-l):vswf.jlm(l,l)]*
            eimphi
         l<-l-1
#-------------------------------------------------------------------------------
# SSH DESLOCADOS
#-------------------------------------------------------------------------------
         Y.l.m[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l-1,mm-1)]
         Y.l.o[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l-1,mm  )]
         Y.l.M[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l-1,mm+1)]
#------------------------------------------------
         Y.o.m[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l  ,mm-1)]
         Y.o.o[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l  ,mm  )]
         Y.o.M[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l  ,mm+1)]
#------------------------------------------------
         Y.L.m[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l+1,mm-1)]
         Y.L.o[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l+1,mm  )]
         Y.L.M[vswf.jlm(l,mm)]<-Ylm[vswf.jlm(l+1,mm+1)]
#------------------------------------------------
         c.l[vswf.jlm(l,mm)]<-cl(l,mm)
         l.l[vswf.jlm(l,mm)]<-ll(l,mm)
#------------------------------------------------
         c.p[vswf.jlm(l,mm)]<-cp(l,mm)
         c.m[vswf.jlm(l,mm)]<-cm(l,mm)
         c.z[vswf.jlm(l,mm)]<-cz(l,mm)
#------------------------------------------------
         KmLm[vswf.jlm(l,mm)]<-Km(l+1,mm,mm-1)
         KoLo[vswf.jlm(l,mm)]<-Ko(l+1,mm,mm  )
         KpLM[vswf.jlm(l,mm)]<-Kp(l+1,mm,mm+1)
#------------------------------------------------
         KmlM[vswf.jlm(l,mm)]<-Km(l,mm,mm+1)
         Kolo[vswf.jlm(l,mm)]<-Ko(l,mm,mm  )
         Kplm[vswf.jlm(l,mm)]<-Kp(l,mm,mm-1)
#------------------------------------------------
         Km.lm[vswf.jlm(l,mm)]<-KYmlm(l,mm)
         Ko.lm[vswf.jlm(l,mm)]<-KYolm(l,mm)
         Kp.lm[vswf.jlm(l,mm)]<-KYplm(l,mm)
#------------------------------------------------
         Km.Lm[vswf.jlm(l,mm)]<-KYmLm(l,mm)
         Ko.Lm[vswf.jlm(l,mm)]<-KYoLm(l,mm)
         Kp.Lm[vswf.jlm(l,mm)]<-KYpLm(l,mm)
#------------------------------------------------
      }
#-------------------------------------------------------------------------------
# VECTOR SPHERICAL HARMONICS
#-------------------------------------------------------------------------------
      # R VERSOR (\bm{{\rm\hat{r}}})
      rm<-rep(sth*(cph-1i*sph)/sqrt(2),LMAX)
      rz<-rep(cth,LMAX)
      rp<-rep(sth*(cph+1i*sph)/sqrt(2),LMAX)
#------------------------------------------------
      # X - EXPLICIT EXPRESSION
      Xm<-l.l*c.m*Y.o.m/sqrt(2)
      Xz<-l.l*c.z*Y.o.o
      Xp<-l.l*c.p*Y.o.M/sqrt(2)
#------------------------------------------------
      # Y - EXPLICIT EXPRESSION
      Ym<-(-Kplm*Y.l.m+KmLm*Y.L.m)/sqrt(2)
      Yz<-  Kolo*Y.l.o+KoLo*Y.L.o
      Yp<-( KmlM*Y.l.M-KpLM*Y.L.M)/sqrt(2)
#------------------------------------------------
      # V - EXPLICIT EXPRESSION
      Vm<-l.l*(-(c.l+1)*Kplm*Y.l.m-c.l*KmLm*Y.L.m)/sqrt(2)
      Vz<-l.l*( (c.l+1)*Kolo*Y.l.o-c.l*KoLo*Y.L.o)
      Vp<-l.l*( (c.l+1)*KmlM*Y.l.M+c.l*KpLM*Y.L.M)/sqrt(2)
#------------------------------------------------
      # Yl,l-1,m  (BLATT AND WEISSKOPF)
      Ym.lm<--Km.lm*Y.l.m
      Yz.lm<- Ko.lm*Y.l.o
      Yp.lm<- Kp.lm*Y.l.M
#------------------------------------------------
      # Yl,l+1,m  (BLATT AND WEISSKOPF)
      Ym.Lm<--Km.Lm*Y.L.m
      Yz.Lm<--Ko.Lm*Y.L.o
      Yp.Lm<- Kp.Lm*Y.L.M 
##------------------------------------------------
   #   # Y - DEFINED BY r Y_lm (DIRECT MULTIPLICATION)
   #   Ym.o<-rm*Ylm[1:LMAX]
   #   Yz.o<-rz*Ylm[1:LMAX]
   #   Yp.o<-rp*Ylm[1:LMAX]
   #------------------------------------------------
   #   # V - DEFINED BY -ir x X (CROSS PRODUCT)
   #   Vm.o<-rm*Xz-rz*Xm
   #   Vz.o<-rp*Xm-rm*Xp
   #   Vp.o<-rz*Xp-rp*Xz
   #------------------------------------------------
   #   # Y - DEFINED BY BLATT & AND WEISSKOPF
   #   Ym.y<-(sqrt(c.l)*Ym.lm-sqrt(c.l+1)*Ym.Lm)/sqrt(2*c.l+1)
   #   Yz.y<-(sqrt(c.l)*Yz.lm-sqrt(c.l+1)*Yz.Lm)/sqrt(2*c.l+1)
   #   Yp.y<-(sqrt(c.l)*Yp.lm-sqrt(c.l+1)*Yp.Lm)/sqrt(2*c.l+1)
   ##------------------------------------------------
   #   # V - DEFINED BY BLATT & WEISSKOPF
   #   Vm.y<-(sqrt(c.l+1)*Ym.lm+sqrt(c.l)*Ym.Lm)/sqrt(2*c.l+1)
   #   Vz.y<-(sqrt(c.l+1)*Yz.lm+sqrt(c.l)*Yz.Lm)/sqrt(2*c.l+1)
   #   Vp.y<-(sqrt(c.l+1)*Yp.lm+sqrt(c.l)*Yp.Lm)/sqrt(2*c.l+1)
##-------------------------------------------------------------------------------
   #   # Y_{l,l-1}^m - BLATT & WEISSKOPF FROM XVY VSH
   #   Ym.lm.o<-(sqrt(c.l+1)*Vm+sqrt(c.l)*Ym)/sqrt(2*c.l+1) #OK
   #   Yz.lm.o<-(sqrt(c.l+1)*Vz+sqrt(c.l)*Yz)/sqrt(2*c.l+1) #OK
   #   Yp.lm.o<-(sqrt(c.l+1)*Vp+sqrt(c.l)*Yp)/sqrt(2*c.l+1) #OK
##-------------------------------------------------------------------------------
   #   # Y_{l,l+1}^m - BLATT & WEISSKOPF FROM XVY VSH
   #   Ym.Lm.o<-(sqrt(c.l)*Vm-sqrt(c.l+1)*Ym)/sqrt(2*c.l+1) #OK
   #   Yz.Lm.o<-(sqrt(c.l)*Vz-sqrt(c.l+1)*Yz)/sqrt(2*c.l+1) #OK
   #   Yp.Lm.o<-(sqrt(c.l)*Vp-sqrt(c.l+1)*Yp)/sqrt(2*c.l+1) #OK
#-------------------------------------------------------------------------------
# RESULTADOS
#-------------------------------------------------------------------------------
      Ylm<-Ylm[1:LMAX]
      Ym[1]<-Ylm[1]*sth*(cph-1i*sph)/sqrt(2)
      Yz[1]<-Ylm[1]*cth
      Yp[1]<-Ylm[1]*sth*(cph+1i*sph)/sqrt(2)
      U<-data.frame(l=c.l,m=c.z,
   #                  r=l.l,cm=c.m,cp=c.p,
   #                  KmLm,KoLo,KpLM,KmlM,Kolo,Kplm,
                     Y=Ylm, #SCALAR SPHERICAL HARMONICS
                     VM=Vm,VZ=Vz,VP=Vp, # V VECTOR SPHERICAL HARMONIC
                     XM=Xm,XZ=Xz,XP=Xp, # X VECTOR SPHERICAL HARMONIC
                     YM=Ym,YZ=Yz,YP=Yp  # Y VECTOR SPHERICAL HARMONIC
##                     Ym.lm,Yz.lm,Yp.lm, # B&W Y_{l,l-1}^m
##                     Ym.Lm,Yz.Lm,Yp.Lm  # B&W Y_{l,l+1}^m
                     #Ym.y,Yz.y,Yp.y,
                     #Vm.y,Vz.y,Vp.y,
                     #Ym.o,Yz.o,Yp.o,
                     #Vm.o,Vz.o,Vp.o,
                     #Ym.lm.o,Yz.lm.o,Yp.lm.o,
                     #Ym.Lm.o,Yz.Lm.o,Yp.Lm.o
      )
      return(U)
   }
}
