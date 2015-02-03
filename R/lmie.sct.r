#' Scattering coeficients in the Lorentz-Mie
#' 
#' @details Scattering and Extinction coefficients.
#' @param m Ratio of refraction indices.
#' @param x Form factor
#' @param dist If is to distinguish contributions of \eqn{a_n} and \eqn{b_n}.
#' @return Extinction and Scattering coefficients \eqn{Qsca} and \eqn{Qext} if
#' \code{dist=FALSE}. If \code{dist=TRUE}, so the contribution of \eqn{a_n} and
#' \eqn{b_n} are distinguished in \eqn{Q_{sca}\to Q_{sca}+Q_{scb}} and 
#' \eqn{Q_{ext}\to Q_{exa}+Q_{exb}}.
#' @import lmie.exp
#' @export
#' @examples
#' dx<-.1   # FAST VALUE
#' #dx<-.001 # HIGH RESOLUTION
#' m<-1.33
#' x<-seq(dx,50,dx)
#' u<-lmie.sct(m,x)
#' v<-lmie.sct(m,x,dist=TRUE)
#' plot(x,u$Qsca,type='l')
#' points(x,v$Qsca,type='l',col='red')
#' points(x,v$Qscb,type='l',col='blue')
lmie.sct<-function(m,x,dist=FALSE){
   if(length(x)>1){
      Q<-sapply(x,lmie.sct,m=m,dist=dist)
      Q<-as.data.frame(apply(t(Q),2,as.numeric))
      
      return(Q)
   }
   ab<-lmie.exp(m,x)
   n<-nrow(ab)
   n<-1:n
   an<-ab$an
   bn<-ab$bn
   if(dist){
      Qsca<-(2/(x^2))*sum((2*n+1)*abs(an)^2)
      Qscb<-(2/(x^2))*sum((2*n+1)*abs(bn)^2)
      Qexa<-(2/(x^2))*sum((2*n+1)*Re(an))
      Qexb<-(2/(x^2))*sum((2*n+1)*Re(bn))
      Q<-data.frame(Qsca,Qscb,Qexa,Qexb)
   }else{
      Qsca<-(2/(x^2))*sum((2*n+1)*(abs(an)^2+abs(bn)^2))
      Qext<-(2/(x^2))*sum((2*n+1)*Re(an+bn))
      Q<-data.frame(Qsca,Qext)
   }
   return(Q)
}

