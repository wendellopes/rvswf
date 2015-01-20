#' Auxiliary function used in calculations of Bessel Beams and Cylindrical Wave Guides.
#' 
#' @details The definition is given by 
#' \eqn{\psi_m(\bm{k};\bm{r})=J_m(kr)e^{ism\phi}e^{ik_zz}}.
#' @param gama The transversal component of the wave vector.
#' @param kz   The longitudinal component of the wave vector.
#' @param  x The component \eqn{x} of the position vector.
#' @param  y The component \eqn{y} of the position vector.
#' @param  z The component \eqn{z} of the position vector.
#' @param  s The chirality of the function.
#' @return The complex value of \eqn{\psi_m(x)}.
#' @examples
#' vswf.psi(1,2,3,4,5,6)
vswf.psi<-function(gama,kz,x,y,z,m,s=1){
   rho<-sqrt(x^2+y^2)
   if(abs(rho)<1e-14){
      cph<-1
      sph<-0
   }else{
      cph<-x/rho
      sph<-y/rho   
   }
   u<-besselJ(gama*rho,m)*((cph+1i*s*sph)^m)*exp(1i*kz*z)
   return(u)
}
