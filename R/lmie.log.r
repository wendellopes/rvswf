#' Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.
#' 
#' @details Calculation done using the logarithmic derivative of
#' Ricatti-Bessel functions.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param NMAX The maximum value of \code{n}. Default value from paper.
#' @return The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and also
#' the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}. The maximum value of
#' \code{n} wil
#' @examples
#' m<-1.2+.1i
#' x<-3
#' print(lmie.rho(1.2,5))
#' 
#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF LOGARITHMIC DERIVATIVES
#-------------------------------------------------------------------------------
lmie.log<-function(m,x,NMAX=floor(abs(x+7.5*x^(1/3))+2)){
   # CALCULATIONS - It must exclude first value (n=0)
   An.1<-lcfa.ric(NMAX  ,x)$Cn[-1]
   An.m<-lcfa.ric(NMAX,m*x)$Cn[-1]
   #------------------------------------
   # UPWARD RECURRENCE
   Cn<-Bn<-rep(1,NMAX)
   Cn[1]<-1/(1+1i*(cos(x)+x*sin(x))/(sin(x)-x*cos(x)))
   Bn[1]=-lcfe.afs(1,x)+1/(lcfe.afs(1,x)+1i)
   for(n in 2:NMAX){
      Bn[n]<--lcfe.afs(n,x)+1/(lcfe.afs(n,x)-Bn[n-1])
      Cn[n]<-Cn[n-1]*(Bn[n]+n/x)/(An.1[n]+n/x)
   }
   # OTHER ExPRESSIONS
   Ta<-(An.m/m-An.1)/(An.m/m-Bn)
   Tb<-(An.m*m-An.1)/(An.m*m-Bn)
   # MIE COEFFICIENTS
   an<-Cn*Ta
   bn<-Cn*Tb
   # SCATTERING
   u<-data.frame(Cn,Ta,Tb,an,bn)
   return(u)
}
