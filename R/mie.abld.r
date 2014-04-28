#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF LOGARITHMIC DERIVATIVES
#-------------------------------------------------------------------------------
mie.abld<-function(m,x,NMAX=floor(x+7.5*x^(1/3))+2){
   # STOPPING CRITERIUM
   # DOWNWARD RECURRENCE
   An.1<-rep(-17,NMAX)
   An.m<-rep(-17,NMAX)
   # CALCULATIONS
   An.1[NMAX]<-lcf.rbld(NMAX,  x)
   An.m[NMAX]<-lcf.rbld(NMAX,m*x)
   for(n in NMAX:2){
      An.1[n-1]<-lcf.afsn(n,  x)-1/(n/(  x)+An.1[n])
      An.m[n-1]<-lcf.afsn(n,m*x)-1/(n/(m*x)+An.m[n])
   }
   # UPWARD RECURRENCE
   Cn<-Bn<-rep(1,NMAX)
   Cn[1]<-1/(1+1i*(cos(x)+x*sin(x))/(sin(x)-x*cos(x)))
   Bn[1]=-lcf.afsn(1,x)+1/(lcf.afsn(1,x)+1i)
   for(n in 2:NMAX){
      Bn[n]<--lcf.afsn(n,x)+1/(lcf.afsn(n,x)-Bn[n-1])
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
