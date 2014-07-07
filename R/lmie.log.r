#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF LOGARITHMIC DERIVATIVES
#-------------------------------------------------------------------------------
lmie.log<-function(m,x,NMAX=floor(x+7.5*x^(1/3))+2){
   # STOPPING CRITERIUM
   # DOWNWARD RECURRENCE
   An.1<-rep(-17,NMAX)
   An.m<-rep(-17,NMAX)
   #------------------------------------
   # CALCULATIONS - ORIG
   #An.1[NMAX]<-lcfe.rbl(NMAX,  x)
   #An.m[NMAX]<-lcfe.rbl(NMAX,m*x)
   #for(n in NMAX:2){
   #   An.1[n-1]<-lcfe.afs(n,  x)-1/(n/(  x)+An.1[n])
   #   An.m[n-1]<-lcfe.afs(n,m*x)-1/(n/(m*x)+An.m[n])
   #}
   #------------------------------------
   # CALCULATIONS - NEW
   An.1<-lcfa.ric(NMAX  ,x)$Dn
   An.1<-lcfa.ric(NMAX,m*x)$Dn
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
