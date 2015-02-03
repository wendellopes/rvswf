#' Beam Shape Coefficients for a Generic Plane Wave
#' 
#' @details Calculate the Beam Shape Coefficients for a plane wave with
#' arbitrary polarization and propagating in arbitrary direction. It can
#' be compared with the specific case calculated by Mie or that done in
#' the book of Jackson, called here by Mie Plane Wave.
#' @param kx The component \eqn{x} of the wave vector.
#' @param ky The component \eqn{y} of the wave vector.
#' @param kz The component \eqn{z} of the wave vector.
#' @param  x The component \eqn{x} of the position vector.
#' @param  y The component \eqn{y} of the position vector.
#' @param  z The component \eqn{z} of the position vector.
#' @param xo The component \eqn{x} of the origin vector.
#' @param yo The component \eqn{y} of the origin vector.
#' @param zo The component \eqn{z} of the origin vector.
#' @param implicit The way it calculates the Beam Shape Coefficients. If TRUE,
#' so we use the function \code{vswf.vsh} to calculate. If FALSE, we use
#' the explicit expression for the components.
#' @include vswf.vsh.r vswf.qlm.r vswf.jlm.r
#' @export
#' @seealso \code{\link{vswf.vsh}}, \code{\link{vswf.mpw}}, \code{\link{vswf.qlm}},
#' \code{\link{vswf.jlm}}.
#' @examples
#' lm<-5
#' a<-vswf.mpw(lm,norm=TRUE)
#' b<-vswf.gpw(0,0,1,1/sqrt(2),1i/sqrt(2),0,lm)
#' #
#' plot(Re(a$GTE),type='b')
#' points(Re(b$GTE),pch=4,col='red',type='b')
#' #
#' plot(Re(a$GTM),type='b')
#' points(Re(b$GTM),pch=4,col='red',type='b')
#' #
#' plot(Im(a$GTE),type='b')
#' points(Im(b$GTE),pch=4,col='red',type='b')
#' #
#' plot(Im(a$GTM),type='b')
#' points(Im(b$GTM),pch=4,col='red',type='b')
vswf.gpw<-function(kx,ky,kz,ux,uy,uz,lmax,xo=0,yo=0,zo=0,implicit=TRUE){
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
   # VECTORS
   GTE<-rep(0,LMAX)
   GTM<-rep(0,LMAX)
   Il<-rep(0,LMAX)
   # Implicit calculation
   if(implicit){
      for(l in 1:lmax){
         m<--l:l
         lmz<-vswf.jlm(l,m)
         Il[lmz]<-1i^l
      }
      GTE<-4*pi*Il*(Conj(a$Xm)*hum+Conj(a$Xz)*huz+Conj(a$Xp)*hup)
      GTM<-4*pi*Il*(Conj(a$Xm)*ham+Conj(a$Xz)*haz+Conj(a$Xp)*hap)
   }else{
      U<-vswf.qlm(kz/k,lmax+1)
      qlm<-U$Qlm
      rm(U)
      # Constants for Normalized Associated Legendre Polynomials
      cplm<-function(l,m){return(sqrt(l*(l+1)-m*(m+1)))}
      cmlm<-function(l,m){return(sqrt(l*(l+1)-m*(m-1)))}
      il<-function(l){return((1i^l)/sqrt(l*(l+1)))}
      # Calculations
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
         GTE[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hup/sqrt(2)+
            m*ez*qlm[lmz]*huz+cmlm(l,m)*em*qlm[lmm]*hum/sqrt(2))
         GTM[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hap/sqrt(2)+
            m*ez*qlm[lmz]*haz+cmlm(l,m)*em*qlm[lmm]*ham/sqrt(2))
      }
   } 
   return(data.frame(GTE,GTM))
}
