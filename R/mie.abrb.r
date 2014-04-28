#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF RATIO BETWEEN BESSEL FUNCTIONS
#-------------------------------------------------------------------------------
mie.abrb<-function(m,x,NMAX=floor(x+7.5*x^(1/3))+2,DIRECT=TRUE){
   # STOPPING CRITERIUM
   # We need the zeroth therm at postion [1]
   rho.1<-rho.m<-g<-Cn<-rep(-17,NMAX+1)
   # Calculatin \gamma_0=g[1]
      p0<-sin(x)
      q0<-cos(x)
      z0<-p0+1i*q0
      p1<-p0/x-q0
      q1<-q0/x+p0
      z1<-p1+1i*q1
   # Starting series
   g[1]<-z0/z1
   Cn[1]<-p0/z0
   # STARTING VALUES 
   if(DIRECT){
      rho.1[NMAX+1]<-lcf.sbrd(NMAX,  x)
      rho.m[NMAX+1]<-lcf.sbrd(NMAX,m*x)
   }else{
      rho.1[NMAX+1]<-1/lcf.sbri(NMAX,  x)
      rho.m[NMAX+1]<-1/lcf.sbri(NMAX,m*x)
   }
   # DOWNWARD RECURRENCE
   for(n in NMAX:1){
      rho.1[n]<-lcf.afsn(2*n+1,  x)-1/rho.1[n+1]
      rho.m[n]<-lcf.afsn(2*n+1,m*x)-1/rho.m[n+1]
   }
   # UPWARD RECURRENE 
   for(n in 1:NMAX){
      g[n+1]<-1/(lcf.afsn(2*n+1,x)-g[n])
      Cn[n+1]<-Cn[n]*g[n]/rho.1[n]
   }
   # OTHER ExPRESSIONS
   n<-1:(NMAX+1)
   k<-(1-1/m^2)*n/x
   Ta<-(rho.m/m-rho.1+k)/(rho.m/m-g+k)
   Tb<-(rho.m*m-rho.1  )/(rho.m*m-g  )
   # MIE COEFFICIENTS
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
