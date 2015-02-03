#' Beam Shape Coefficient for Transversal Bessel Beam.
#' 
#' @details Calculate the Beam Shape Coefficients for the partial wave expansion
#' in Vector Spherical Wave Functions.
#' @param TE Type of the wave field (TE or TM).
#' @param gama The transversal component of the wave vector \eqn{\gamma}.
#' @param  kz  The longitudinal component of the wave vector \eqn{k_z}.
#' @param xo The position \eqn{x} at which the expansion is performed.
#' @param yo The position \eqn{y} at which the expansion is performed.
#' @param zo The position \eqn{z} at which the expansion is performed.
#' @param lmax The maximum value of \eqn{l} in the expansion.
#' @param M Order of the Bessel Beam.
#' @param s Chirality of the Bessel Beam (\eqn{S=\pm 1}).
#' @seealso \code{\link{vwfd.bbp}}, \code{\link{vswf.bbp}}, \code{\link{vswf.pwe}}.
#' @include vswf.qlm, vswf.jlm, vswf.psi
#' @export
#' @examples
#' lambda=.5e-6            # Propagating wavelength
#' a<-8*lambda             # Size x of the waveguide (Rectangular Wave Guide)
#' b<-8*lambda             # Size y of the waveguide (Rectangular Wave Guide)
#' l<-6*lambda             # Size z of the waveguide (Rectangular Wave Guide)
#' M<-3                    # x wavefield mode
#' N<-5                    # y wavefield mode
#' S<--1                   # Chirality
#' TM<-FALSE               # Mode (TM-TRUE/TE-FALSE)
#' # Wave Field Parameters
#' k=2*pi/lambda           # Propagating wavenumber
#' kx<-M*pi/a              # x component of the wavevector
#' ky<-N*pi/b              # y component of the wavevector
#' gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
#' kz<-sqrt(k^2-gama^2)    # z component of the wavevector
#' #-------------------------------------------------------------------------------
#' # Geometry of the calculations
#' #-------------------------------------------------------------------------------
#' NPX<-4*50+1
#' NPY<-4*60+1
#' NPZ<-4*20+1
#' dx<-2*a/(NPX-1)
#' dy<-2*b/(NPY-1)
#' dz<-2*l/(NPZ-1)
#' x<-seq(-a,a,by=dx)
#' y<-seq(-b,b,by=dy)
#' z<-0
#' #-------------------------------------------------------------------------------
#' lmax<-80 # Number of partial waves to sum (for some l we have 2l+1 values of m)
#' #-------------------------------------------------------------------------------
#' # POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#' #-------------------------------------------------------------------------------
#' xo<-yo<-zo<-0
#' #-------------------------------------------------------------------------------
#' # CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#' #-------------------------------------------------------------------------------
#' x<-sample(x,1)
#' y<-sample(y,1)
#' z<-0
#' #-------------------------------------------------------------------------------
#' # BSC CALCULATIONS
#' #-------------------------------------------------------------------------------
#' BBZ<-vwfd.bbz(TM,M,S,gama,kz,x+xo,y+yo,z+zo)
#' BSC<-vswf.bbz(TM,gama,kz,xo,yo,zo,lmax,M,S)
#' PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
#' #-------------------------------------------------------------------------------
#' # VALUES
#' #-------------------------------------------------------------------------------
#' cat("Distance x from origin in wavelength (from ",-a/lambda,"to ",a/lambda,"):",x/lambda,"\n")
#' cat("Distance y from origin in wavelength (from ",-b/lambda,"to ",b/lambda,"):",y/lambda,"\n")
#' df<-data.frame(
#'    PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
#'    BBZ=c(BBZ$Em,BBZ$Ez,BBZ$Ep,BBZ$Hm,BBZ$Hz,BBZ$Hp),
#'    row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
#'    )
#' df$DIF<-df$PWE-df$BBZ
#' print(df)
vswf.bbz<-function(TM=TRUE,gama,kz,xo,yo,zo,lmax,M,s){
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
   if(TM){
      return(data.frame(GTE,GTM))
   }else{
      return(data.frame(GTE=GTM,GTM=-GTE))
   }
}
