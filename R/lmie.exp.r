#' Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.
#' 
#' @details Parser to the two ways to calculate the Lorentz-Mie coefficents.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param by \code{LD} using logarithmic derivative; \code{RB} using the ratio
#' between Ricatti-Bessel functions.
#' @return The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and also
#' the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}.
#' @examples
#' m<-5
#' x<-3
#' a<-lmie.exp(m,x,by="LD")
#' b<-lmie.exp(m,x,by="RB")
#' print(a-b) # Must be close to zero
lmie.exp<-function(m,x,by="LD",...){
   if(!by%in%c("LD","RB")){
      stop("by must be Logorithm Derivative (LD) or 
             Ratio between Bessel functions (RB)")
   }
   if(by=="LD"){
   	return(lmie.log(m,x,...))
   }else{
      return(lmie.rho(m,x,...))
   }
}
