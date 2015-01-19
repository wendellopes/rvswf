#' Convert the value of (l,m) to a single value
#' 
#' @param l The orbital angular momentum eigenvalue
#' @param m The eigenvalue of the projection of angular momentum on \eqn{z} axis
#' @return A unique value used to transform a double summation on a single one.
#' @examples 
#' for(l in 0:5){
#'    for(m in -l:l){
#'       print(c(l,m,vswf.jlm(l,m)))
#'    }
#' }
vswf.jlm<-function(l,m){
   j.lm<-l*(l+1)+m+1    
   j.tf<-abs(m)>l
   j.lm[j.tf]<-1
   return(j.lm)
}
