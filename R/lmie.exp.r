#' Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.
#' 
#' @details Parser to the two ways to calculate the Lorentz-Mie coefficents.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param by \code{LD} using logarithmic derivative; \code{RB} using the ratio
#' between Ricatti-Bessel functions.
#' @return The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and also
#' the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}.
#' @include lmie.log, lmie.rho
#' @export
#' @examples
#' # Table 4.1 - Bohren and Hoffman Book, pg 114
#' # The values showed bellow are conjugated. Probably
#' # there is tipography error in the book.
#' m<-1.33+1e-8i
#' x<-3
#' n<-c(1,2,3,4,5)
#' an<-c(5.1631e-1+4.9973e-1i, 3.4192e-1+4.7435e-1i,4.8467e-2+2.1475e-1i,
#'       1.0346e-3+3.2148e-2i, 9.0375e-6+3.0062e-3i)
#' bn<-c(7.3767e-1+4.3990e-1i, 4.0079e-1+4.9006e-1i,9.3553e-3+9.6269e-2i,
#'       6.8810e-5+8.2949e-3i, 2.8390e-7+5.3204e-4i)
#' bh<-data.frame(n,an,bn)
#' a<-lmie.exp(m,x,by="LD",NMAX=5)
#' b<-lmie.exp(m,x,by="RB",NMAX=5)
#' print(a)
#' print(b)
#' print(bh)
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
