#' Beam Shape Coefficient for Circularly polarized Bessel Beam
#' 
#' @details Calculate the Beam Shape Coefficients for the partial wave expansion
#' in Vector Spherical Wave Functions.
#' @param gama The transversal component of the wave vector \eqn{\gamma}.
#' @param  kz  The longitudinal component of the wave vector \eqn{k_z}.
#' @param xo The position \eqn{x} at which the expansion is performed.
#' @param yo The position \eqn{y} at which the expansion is performed.
#' @param zo The position \eqn{z} at which the expansion is performed.
#' @param lmax The maximum value of \eqn{l} in the expansion.
#' @param M Order of the Bessel Beam.
#' @param s Chirality of the Bessel Beam (\eqn{S=\pm 1}).
#' @param p Polarity of the Bessel beam (\eqn{S=\pm 1}).
#' @return The Beam Shape Coefficients for a Circularly polarized Bessel Beam.
#' @seealso \code{\link{vwfd.bbp}}, \code{\link{vswf.bbp}}, \code{\link{vswf.pwd}}.
#' @export
#' @examples
#' lambda=.5e-6            # Propagating wavelength
#' a<-8*lambda             # Size x of the waveguide (Rectangular Wave Guide)
#' b<-8*lambda             # Size y of the waveguide (Rectangular Wave Guide)
#' l<-6*lambda             # Size z of the waveguide (Rectangular Wave Guide)
#' M<-5                    # x wavefield mode
#' N<-4                    # y wavefield mode
#' S<--1                   # Chirality
#' P<- 1                    # Mode (+/- 1)
#' # Wave Field Parameters
#' k=2*pi/lambda           # Propagating wavenumber
#' kx<-M*pi/a              # x component of the wavevector
#' ky<-N*pi/b              # y component of the wavevector
#' gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
#' kz<-sqrt(k^2-gama^2)    # z component of the wavevector
#' TM<-ifelse(P==-1,TRUE,FALSE)
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
#' lmax<- 80 # Number of partial waves to sum (for some l we have 2l+1 values of m)
#' #-------------------------------------------------------------------------------
#' # POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#' #-------------------------------------------------------------------------------
#' # ARBITRARY: We select two points, the is for the origin of the expansion,
#' #            the second is for the test, that must have be move to the new
#' #            system of coordinates x'=x-xo
#' xo<-yo<-zo<-0
#' x<-sample(x,1)
#' y<-sample(y,1)
#' #-------------------------------------------------------------------------------
#' # BSC CALCULATIONS
#' #-------------------------------------------------------------------------------
#' BBP<-vwfd.bbp(P,M,S,gama,kz,x+xo,y+yo,z+zo)  # WAVE FIELD DEFINITION
#' BSC<-vswf.bbp(gama,kz,xo,yo,zo,lmax,M,S,P)   # BEAM SHAPE COEFFICIENTS
#' PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)  # PARTIAL WAVE EXPANSION
#' #-------------------------------------------------------------------------------
#' # VALUES
#' #-------------------------------------------------------------------------------
#' cat("Distance x from origin in wavelength (from ",-a/lambda,"to ",a/lambda,"):",x/lambda,"\n")
#' cat("Distance y from origin in wavelength (from ",-b/lambda,"to ",b/lambda,"):",y/lambda,"\n")
#' df<-data.frame(
#'    PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
#'    BBP=c(BBP$Em,BBP$Ez,BBP$Ep,BBP$Hm,BBP$Hz,BBP$Hp),
#'    row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
#'    )
#' df$DIF<-df$PWE-df$BBP
#' print(df)
vswf.bbp<-function(gama,kz,xo,yo,zo,lmax,M,s,p,code="C"){
   LMAX=lmax*(lmax+2)+1
   dummy<-rep(0,LMAX)
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      u<-.C("vswf_bbp",
         M=as.integer(M),
         s=as.integer(s),
         p=as.integer(p),
         lmax=as.integer(lmax),
         NMAX=as.integer(200),

         gamma=as.double(gama),
         kz=as.double(kz),

         x=as.double(x),
         y=as.double(y),
         z=as.double(z),

         GTE=as.complex(dummy),      # an
         GTM=as.complex(dummy))      # bn
      return(data.frame(GTE=u$GTE,GTM=u$GTM))
   }else{
      k<-sqrt(gama^2+kz^2)
      cth<-kz/k
      sth<-gama/k
      #----------------------------------------
      u<-vswf.qlm(kz/k,lmax+1)
      Qlm<-u$Qlm
      dQlm<-u$dQlm
      rm(u)
      #----------------------------------------
      clmp<-function(l,m,p){return(sqrt(l*(l+1)-m*(m+p)))}
      Klmqp<-function(l,m,q,p){return(((l+p*m)*(l+p*q))/((2*l-1)*(2*l+1)))}
      #----------------------------------------
      GTE<-rep(0,LMAX)
      GTM<-rep(0,LMAX)
      for(l in 1:lmax){
         m<--l:l
         psio<-vswf.psi(gama,kz,xo,yo,zo,M-s*(m-p),s)
         if(l>0){
            A0<-4*pi*(1i^(l-s*(m-p)))*psio/sqrt(2*l*(l+1))
         }else{
            A0<-0
         }
         GTE[vswf.jlm(l,m)]<-A0*clmp(l,m,-p)*Qlm[vswf.jlm(l,m-p)]
         GTM[vswf.jlm(l,m)]<--1i*A0*(cth*sth*dQlm[vswf.jlm(l,m)]-
            p*m*Qlm[vswf.jlm(l,m)]/sth)  
      }
      return(data.frame(GTE,GTM))
   }
}
