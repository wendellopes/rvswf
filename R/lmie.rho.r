#' Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.
#' 
#' @details Calculation done using the ratio between Ricatti-Bessel functions.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param NMAX The maximum value of \code{n}. Default value from paper.
#' @return The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and also
#' the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}.
#' @export
#' @examples
#' n<-5
#' x<-3
#' print(lmie.rho(1.2,5))
#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF RATIO BETWEEN BESSEL FUNCTIONS
#-------------------------------------------------------------------------------
lmie.rho<-function(m,x,NMAX=floor(abs(x+7.5*x^(1/3))+2)){#,DIRECT=TRUE){
   # We need the zeroth therm at postion [1]
   rho.1<-rho.m<-g<-Cn<-rep(-17,NMAX+1)
   p0<-sin(x)
   q0<-cos(x)
   z0<-p0+1i*q0
   p1<-p0/x-q0
   q1<-q0/x+p0
   z1<-p1+1i*q1
   # Starting series
   g[1]<-z0/z1
   Cn[1]<-p0/z0
   #------------------------------------
   rho.1<-lcfa.ric(NMAX,  x)$rn
   rho.m<-lcfa.ric(NMAX,m*x)$rn
   #------------------------------------
   # UPWARD RECURRENE 
   for(n in 1:NMAX){
      g[n+1]<-1/(lcfe.afs(2*n+1,x)-g[n])
      Cn[n+1]<-Cn[n]*g[n]/rho.1[n]
   }
   # OTHER ExPRESSIONS
   n<-1:(NMAX+1) #n=n+1?
   k<-(1-1/m^2)*n/x
   Ta<-(rho.m/m-rho.1+k)/(rho.m/m-g+k)
   Tb<-(rho.m*m-rho.1  )/(rho.m*m-g  )
   # MIE COEFFICIENTS
   # CALCULATIONS - It must exclude first value (n=0)
   n<-1:NMAX
   Cn<-Cn[n+1]
   Ta<-Ta[n]
   Tb<-Tb[n]
   an<-Cn*Ta
   bn<-Cn*Tb
   # SCATTERING 
   u<-data.frame(Cn,Ta,Tb,an,bn)
   #return(u[2:(NMAX+1),])
   return(u)
}
