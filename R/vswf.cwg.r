#' Beam Shape Coefficients for Cylindrical Wave Guides.
#' 
#' @details  Calculates the Beam Shape Coefficients used to do the Partial Wave 
#' Expansion on Vector Spherical Wave Functions.
#' @param TM Type of the wave field.
#' @param gama Transversal component of the wave vector \eqn{\gamma} (single value).
#' @param kz Component \eqn{z} of the wave vector (single value).
#' @param  x Component \eqn{x} of the origin of the expansion (vector).
#' @param  y Component \eqn{y} of the origin of the expansion (vector).
#' @param  z Component \eqn{z} of the origin of the expansion (vector).
#' @param  m Order of the wave guide.
#' @param  s Chirality of the wave guide (\eqn{S=\pm 1})
#' @return The Beam Shape Coefficients \eqn{G^{TE}_{lm}} and \eqn{G^{TM}_{lm}}.
#' @seealso \code{\link{vwfd.cwg}}, \code{\link{vswf.pwe}}, \code{\link{vswf.gwg}}.
#' @examples
#' #-------------------------------------------------------------------------------
#' # Zeros of Bessel Functions
#' #-------------------------------------------------------------------------------
#' #Nth zero (lines) of the Mth Bessel functions (column)
#' ZJ <-matrix(c( 2.4048, 3.8317, 5.1356, 6.3802, 7.5883, 8.7715,
#'                5.5201, 7.0156, 8.4172, 9.7610,11.0647,12.3386,
#'                8.6537,10.1735,11.6198,13.0152,14.3725,15.7002,
#'                11.7915,13.3237,14.7960,16.2235,17.6160,18.9801,
#'                14.9309,16.4706,17.9598,19.4094,20.8269,22.2178),
#'             ncol=6,byrow=TRUE)
#' #-------------------------------------------------------------------------------
#' # Zeros of derivatives of Bessel Functions
#' #-------------------------------------------------------------------------------
#' #Nth zero (lines) of the derivative of Mth Bessel functions (column)
#' ZdJ<-matrix(c( 3.8317, 1.8412, 3.0542, 4.2012, 5.3175, 6.4156,
#'                7.0156, 5.3314, 6.7061, 8.0152, 9.2824,10.5199,
#'                10.1735, 8.5363, 9.9695,11.3459,12.6819,13.9872,
#'                13.3237,11.7060,13.1704,14.5858,15.9641,17.3128,
#'                16.4706,14.8636,16.3475,17.7887,19.1960,20.5755),
#'             ncol=6,byrow=TRUE)
#' #-------------------------------------------------------------------------------
#' # Wave Guide Parameters
#' #-------------------------------------------------------------------------------
#' # Basic Parameters
#' lambda=.5e-6            # Propagating wavelength
#' k=2*pi/lambda           # Propagating wavenumber
#' S=-1                    # Chirality of the beam
#' M<-5  # 0<=M<=5         # Order of the (derivative of the) Bessel function
#' N<-4  # 1<=N<=5         # Order of the zero of the derivative/Bessel function
#' TM<-TRUE                # Mode: TM (Bessel function), TE (derivative)
#' R<-5*lambda             # Radius of the waveguide
#' #-------------------------------------------------------------------------------
#' # MODES: M <- M+1 because J_0(x) ocupies column 1
#' # TE - ZdJ (derivative)
#' # TM - ZJ  (Bessel function)
#' #-------------------------------------------------------------------------------
#' if(TM){
#'    g<-ZJ[N,M+1]/R
#'    MODE<--1
#'    MODE.PWE<-4
#' }else{
#'    g<-ZdJ[N,M+1]/R
#'    MODE<-1
#'    MODE.PWE<-3
#' }
#' kz<-sqrt(k^2-g^2)
#' #-------------------------------------------------------------------------------
#' # Geometry of the calculations
#' #-------------------------------------------------------------------------------
#' NPX<-4*50+1
#' NPY<-4*60+1
#' NPZ<-4*20+1
#' dx<-2*R/(NPX-1)
#' dy<-2*R/(NPY-1)
#' dz<-2*R/(NPZ-1)
#' x<-seq(-R,R,by=dx)
#' y<-seq(-R,R,by=dy)
#' z<-seq(-R,R,by=dz)
#' #-------------------------------------------------------------------------------
#' lmax<- 40  # Number of partial waves to sum (for some l we have 2l+1 values of m)
#' #-------------------------------------------------------------------------------
#' # POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#' #-------------------------------------------------------------------------------
#' xo<-yo<-zo<-0
#' #-------------------------------------------------------------------------------
#' # CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#' #-------------------------------------------------------------------------------
#' x<-x-xo
#' y<-y-yo
#' z<-0
#' x<-sample(x,1)
#' y<-sample(y,1)
#' #-------------------------------------------------------------------------------
#' # BSC CALCULATIONS
#' #-------------------------------------------------------------------------------
#' CWG<-vwfd.cwg(TE=!TM,M,S,g,kz,x+xo,y+yo,z+zo)
#' BSC<-vswf.cwg(TM,g,kz,xo,yo,zo,lmax,m=M,s=S)
#' PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
#' #-------------------------------------------------------------------------------
#' # VALUES
#' #-------------------------------------------------------------------------------
#' cat("Distance x from origin in wavelength (from ",-a/(2*lambda),"to ",a/(2*lambda),"):",x/lambda,"\n")
#' cat("Distance y from origin in wavelength (from ",-b/(2*lambda),"to ",b/(2*lambda),"):",y/lambda,"\n")
#' df<-data.frame(
#'    PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
#'    CWG=c(CWG$Em,CWG$Ez,CWG$Ep,CWG$Hm,CWG$Hz,CWG$Hp),
#'    row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
#'    )
#' df$DIF<-df$PWE-df$CWG
#' print(df)
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
