#' Modified Legendre Polynomials used in scattering theories
#' 
#' @details The functions \eqn{\pi_n(\cos\theta)=\frac{P_n^1(\cos\theta)}{\sin\theta}}
#' and \eqn{\tau_n(\cos\theta)=\frac{d P_n^1(\cos\theta)}{d\theta}}.
#' @param th The angle \eqn{\theta}
#' @param lmax The maximum value of \eqn{l}.
#' @return The functions \eqn{\pi_n} and \eqn{\tau_n}.
#' @export
lmie.leg<-function(th,lmax){
   p<-t<-rep(0,lmax)
   # Initial value
   #p[0]<-0           # \pi_0
   p[1]<-1           # \pi_1
   p[2]<-3*cos(th)   # \pi_2
   #t[0]<-0           # \tau_0
   t[1]<-cos(th)           # \tau_1
   t[2]<-3*cos(2*th) # \tau_2
   for(j in 3:lmax){# j=n+1
      p[j]<-((2*j-1)/(j-1))*cos(th)*p[j-1]-(j/(j-1))*p[j-2]
      t[j]<-j*cos(th)*p[j]-(j+1)*p[j-1]
   }
   return(data.frame(pi=p,tau=t))
}