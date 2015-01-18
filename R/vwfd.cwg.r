#' Vector Wave Field Definition: Cylindrical Wave Guide
#' 
#' @details Calculates the vector components of a cylindrical wave guide of type 
#' Transversal Electric (TE), that means that \eqn{E_z=0}, or Transversal
#' Magnetic (TM), that means \eqn{H_z=0},  in terms of the components of the
#' wave vector and the position vector.
#' @param TE Type of the wave field.
#' @param M Order of the wave guide.
#' @param S Chirality of the wave guide (\eqn{S=\pm 1})
#' @param g Transversal component of the wave vector (\eqn{\gamma}).
#' @param kz Component \eqn{z} of the wave vector (single value).
#' @param x Component \eqn{x} of the position vector (vector).
#' @param y Component \eqn{y} of the position vector (vector).
#' @param z Component \eqn{z} of the position vector (vector).
#' @return A list with the input values, the number of points and the six complex
#' values of components of the electromagnetic fields.
vwfd.cwg<-function(TE=TRUE,M,S,g,kz,x,y,z){
   if(abs(S)!=1){
      stop("S must be plus or minus 1")
   }
if(TE)
{
   te<-1
}else{
   te<-0
}
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("vwfd_cwg",
   TE=as.integer(te),M=as.integer(M),S=as.integer(S),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   gama=as.double(g),kz=as.double(kz),
   x=as.double(x),y=as.double(y),z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
