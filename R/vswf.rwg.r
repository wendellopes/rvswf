#' Beam Shape Coefficients for Rectangular Wave Guides.
#' 
#' @details  Calculates the Beam Shape Coefficients used to do the Partial Wave 
#' Expansion on Vector Spherical Wave Functions.
#' @param TM Type of the wave field.
#' @param kx Component \eqn{x} of the wave vector (single value).
#' @param ky Component \eqn{y} of the wave vector (single value).
#' @param kz Component \eqn{z} of the wave vector (single value).
#' @param  x Component \eqn{x} of the origin of the expansion (vector).
#' @param  y Component \eqn{y} of the origin of the expansion (vector).
#' @param  z Component \eqn{z} of the origin of the expansion (vector).
#' @return The Beam Shape Coefficients \eqn{G^{TE}_{lm}} and \eqn{G^{TM}_{lm}}.
#' @seealso \code{\link{vwfd.rwg}}, \code{\link{vswf.pwe}}, \code{\link{vswf.gwg}}.
#' @examples
#' lambda<-.5e-6           # Propagating wavelength
#' a<-7*lambda             # Size x of the waveguide
#' b<-5*lambda             # Size y of the waveguide
#' M<-6                    # x wavefield mode
#' N<-5                    # y wavefield mode
#' #-------------------------------------------------------------------------------
#' # Wave Field Parameters
#' #-------------------------------------------------------------------------------
#' k<-2*pi/lambda          # Propagating wavenumber
#' kx<-M*pi/a              # x component of the wavevector
#' ky<-N*pi/b              # y component of the wavevector
#' gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
#' kz<-sqrt(k^2-gama^2)    # z component of the wavevector
#' #-------------------------------------------------------------------------------
#' # Geometry of the calculations
#' #-------------------------------------------------------------------------------
#' NPX=200                  # Number of points in each direction (all equal)
#' NPY=200                  # Number of points in each direction (all equal)
#' #-------------------------------------------------------------------------------
#' # Vectors
#' #-------------------------------------------------------------------------------
#' dx<-a/(NPX-1)
#' dy<-b/(NPY-1)
#' x<-seq(0,a,dx)          # x vector of positions
#' y<-seq(0,b,dy)          # y vector of positions
#' z<-0
#' #-------------------------------------------------------------------------------
#' TM<-FALSE
#' lmax<- 40
#' #-------------------------------------------------------------------------------
#' # POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#' #-------------------------------------------------------------------------------
#' # ARBITRARY
#' xo<-a/2
#' yo<-b/2
#' zo<-0
#' #-------------------------------------------------------------------------------
#' # CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#' #-------------------------------------------------------------------------------
#' x<-x-xo
#' y<-y-yo
#' #-------------------------------------------------------------------------------
#' # ARBITRARY POINT FOR CALCULATIONS
#' #-------------------------------------------------------------------------------
#' x<-sample(x,1)
#' y<-sample(y,1)
#' #-------------------------------------------------------------------------------
#' # BSC CALCULATIONS
#' #-------------------------------------------------------------------------------
#' RWG<-vwfd.rwg(TE=!TM,kx,ky,kz,x+xo,y+yo,z+zo)
#' BSC<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,lmax)
#' PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
#' #-------------------------------------------------------------------------------
#' # VALUES
#' #-------------------------------------------------------------------------------
#' cat("Distance x from origin in wavelength (from ",-a/(2*lambda),"to ",a/(2*lambda),"):",x/lambda,"\n")
#' cat("Distance y from origin in wavelength (from ",-b/(2*lambda),"to ",b/(2*lambda),"):",y/lambda,"\n")
#' df<-data.frame(
#'    PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
#'    RWG=c(RWG$Em,RWG$Ez,RWG$Ep,RWG$Hm,RWG$Hz,RWG$Hp),
#'    row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
#'    )
#' df$DIF<-df$PWE-df$RWG
#' print(df)
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
