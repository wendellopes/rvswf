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
#' @include vswf.qlm.r
#' @export
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
#' BSC<-vswf.rwg(kx,ky,kz,xo,yo,zo,lmax,TM)
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
bscf.gpw<-function(M,X,kx,ky,kz,ux,uy,uz,x,y,z,lmax,TM=TRUE){
   x[x==0]<-.Machine$double.xmin
   y[y==0]<-.Machine$double.xmin
   z[z==0]<-.Machine$double.xmin
   nx<-length(x)
   ny<-length(y)
   nz<-length(z)
   dummy<-rep(0,nx*ny*nz)
   tm<-ifelse(TM,1,0)
   u<-.C("bscf_gpw",
      lmax=as.integer(lmax),
      NMAX=as.integer(200),

      kx=as.double(kx),
      ky=as.double(ky),
      kz=as.double(kz),
   
      ux=as.complex(ux),
      uy=as.complex(uy),
      uz=as.complex(uz),
   
      M=as.complex(M),
      X=as.complex(X),

      x=as.double(x),
      y=as.double(y),
      z=as.double(z),

      nx=as.integer(nx),
      ny=as.integer(ny),
      nz=as.integer(nz),

      rx=as.double(dummy),
      ry=as.double(dummy),
      rz=as.double(dummy),

      Fx=as.complex(dummy),
      Fy=as.complex(dummy),
      Fz=as.complex(dummy))
   return(data.frame(x=u$rx,y=u$ry,z=u$rz,u$Fx,u$Fy,u$Fz))
}
