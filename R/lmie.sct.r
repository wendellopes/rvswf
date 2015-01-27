#' Scattering coeficients in the Lorentz-Mie
#' 
#' @details Scattering and Extinction coefficients.
#' @param m Ratio of refraction indices.
#' @param x Form factor
#' @return Extinction and Scattering coefficients
#' @examples
#' # Cantrell papers.
lmie.sct<-function(m,x){
   if(length(x)>1){
      Q<-sapply(x,lmie.sct,m=m)
      Q<-as.data.frame(apply(t(Q),2,as.numeric))
      
      return(Q)
   }
   ab<-lmie.exp(m,x)
   n<-nrow(ab)
   n<-1:n
   an<-ab$an
   bn<-ab$bn
   Qsca<-(2/(x^2))*sum((2*n+1)*(abs(an)^2+abs(bn)^2))
   Qext<-(2/(x^2))*sum((2*n+1)*Re(an+bn))
   Q<-data.frame(Qsca,Qext)
   return(Q)
}

