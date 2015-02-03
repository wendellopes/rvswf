#' Test the calculation
#' 
#' @details Must be rewritten
#' @param n The order of the function
#' @param x The argument of the function
#' @return Data frame with values \eqn{A_n}, \eqn{B_n} and \eqn{C_n}.
#' @include reff.cjn.r reff.cdj.r reff.chn.r reff.cdh.r
#' @export 
test.abc<-function(n,x){
   An<-.5/x+reff.cdj(n+.5,x)/reff.cjn(n+.5,x)
   Bn<-.5/x+reff.cdh(n+.5,x,2)/reff.chn(n+.5,2)
   Cn<-reff.cjn(n+.5,x)/reff.chn(n+.5,x,2)
   return(data.frame(An,Bn,Cn))
}

