#' Vector Wave Field Definition: Circularly polarized Bessel Beam
#' 
#' @details Calculates the vector components of Bessel Beam circularly
#' polarized.
#' @param P Polarity of the Bessel beam (\eqn{S=\pm 1}).
#' @param M Order of the Bessel Beam.
#' @param S Chirality of the Bessel Beam (\eqn{S=\pm 1})
#' @param g Transversal component of the Bessel Beam (\eqn{\gamma}).
#' @param kz Component \eqn{z} of the wave vector (single value).
#' @param x Component \eqn{x} of the position vector (vector).
#' @param y Component \eqn{y} of the position vector (vector).
#' @param z Component \eqn{z} of the position vector (vector).
#' @return A list with the input values, the number of points and the six complex
#' values of components of the electromagnetic fields.
vwfd.bbp<-function(P,M,S,g,kz,x,y,z){
   if(abs(P)!=1){
      stop("P must be plus or minus 1")
   }
   if(abs(S)!=1){
      stop("S must be plus or minus 1")
   }
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("vwfd_bbp",
   P=as.integer(P),M=as.integer(M),S=as.integer(S),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   gama=as.double(g),kz=as.double(kz),
   x=as.double(x),y=as.double(y),z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
