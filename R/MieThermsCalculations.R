#-------------------------------------------------------------------------------
# Calculations must be performed for psi(n,x) and xi(n,x), where
# 1<=n<=N and x=x and x=Mx. Number of points: np=2*N*nx*ny*nz
#-------------------------------------------------------------------------------
# FRACAO CONTINUADA PELO METODO DE LENTZ
# f=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn+...)
# input: vetores an e bn em que an=(a0,a1,...,aN) e bn=(b0,b1,...,bN)
# a0 may be zero or any number
#-------------------------------------------------------------------------------
# BASE NAMES
# lcf : Lentz Continued Fraction Evaluation
# ofc : Optical Force Calculations
# tst : Test Functions
# cmp : Comparison Functions
# aux : Auxiliary Functions
# FUNCTION NAMES
# cbld : Cylindrical Bessel Log Derivative
# cbri : Cylindrical Bessel Ratio Inverse (Barnett)
# cbrd : Cylindrical Bessel Ratio Direct
# sbld : Spherical Bessel Log Derivative
# sbri : Spherical Bessel Ratio Inverse (Barnett)
# sbrd : Spherical Bessel Ratio Direct
# rbld : Riccati-Bessel Log Derivative
# rbri : Riccati-Bessel Ratio Inverse (Barnett)
# rbrd : Riccati-Bessel Ratio Direct
# afsn : Auxiliary function Sn
#-------------------------------------------------------------------------------
# AUXILIARY FUNCTION
lcf.afsn<-function(n,x){
   return(n/x)
}
#-------------------------------------------------------------------------------
# J_{n}/J_{n+1} [OK] DIRECT
lcf.cbrd<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(2*(n+1),x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j+1),x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# j_{n}/j_{n+1} [OK] DIRECT
lcf.sbrd<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(2*(n+1)+1,x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j+1)+1,x); 
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# J_{n+1}/J_{n} [OK] BARNETT
lcf.cbri<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-0    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      an<-(-1)^sign(j)
      j<-j+1
      bn<-lcf.afsn(2*(n+j),x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# j_{n+1}/j_{n} [OK] BARNETT
lcf.sbri<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-0    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      an<-(-1)^sign(j)
      j<-j+1
      bn<-lcf.afsn(2*(n+j)+1,x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# Logarithmic Derivative Cylindrical Bessel [OK]
lcf.cbld<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(n,x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j),x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# Logarithmic Derivative Spherical Bessel [OK]
lcf.sbld<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(n,x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j)+1,x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# Logarithmic Derivative Riccati-Bessel [OK]
lcf.rbld<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(n+1,x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j)+1,x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
#-------------------------------------------------------------------------------
# MIE CODE
# Folowing the notation of Jianqi Shen and Xiaoshu Cai
# Progress in Electromagnetics Research Symposium 2005, Hangzhou, China
# August 22-26
#-------------------------------------------------------------------------------
# MIE CODE
# an=A_n(x)Ta(mx)
# bn=A_n(x)Tb(mx)
# A_n(x)=\psi_n(x)/\xi_n(x)
# Ta(x)=(L_n(mx)/m-L_n(x))/(L_n(mx)/m-B_n(x))
# Tb(x)=(mL_n(mx)-L_n(x))/(mL_n(mx)-B_n(x))
# B_n(x)=\xi_n'(x)/\xi_n(x)
# DOWNWARD RECURRENCE
# L_n(x), L_n(mx)
# L_{n-1}=S_n-1/(S_n+L_n)
# UPWARD RECURRENCE
# A_n(x)=A_{n-1}(x)(B_n(x)+S_n(x))/(L_n(x)+S_n(x))
# B_n(x)=-S_n(x)+1/(S_n(x)-B_{n-1}(x))
# STARTING VALUES
# A_1=1/(1+i(cos(x)+x sin(x))/(sin(x)-x cos(x)))
# B_0=-i
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
#-------------------------------------------------------------------------------
# MIE COEFFICIENTS BY MEANS OF RATIO BETWEEN BESSEL FUNCTIONS
#-------------------------------------------------------------------------------
# TODO: How to calculate \gamma_n=\zeta_{n}/\zeta_{n+1}
#       ** solution: upward recurrence 
# TODO: How to calculate C_n=\psi_n/\zeta_n
#       ** 
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
#-------------------------------------------------------------------------------
mie.anbn<-function(m,x,by="LD",...){
   if(!by%in%c("LD","RB")){
      stop("by must be Logorithm Derivative (LD) or 
             Ratio between Bessel functions (RB)")
   }
   if(by=="LD"){
   	  return(mie.abld(m,x,...))
   }else{
      return(mie.abrb(m,x,...))
   }
}
#-------------------------------------------------------------------------------
# PARTICLE DEPENDENT FORCE CONSTANTS
#-------------------------------------------------------------------------------
ofc.abcn<-function(m,x){
   k<-mie.anbn(m,x)
   n0<-1:(nrow(k)-1)
   n1<-n0+1
   An<-k$an[n1]+Conj(k$an[n0])-2*k$an[n1]*Conj(k$an[n0])
   Bn<-k$bn[n1]+Conj(k$bn[n0])-2*k$bn[n1]*Conj(k$bn[n0])
   Cn<-k$an[n0]+Conj(k$bn[n0])-2*k$an[n1]*Conj(k$an[n0])
   return(data.frame(An,Bn,Cn))
}
#-------------------------------------------------------------------------------
# G^TE, G^TM For Plane Waves +1 (Jackson's Classical Electrodynamics book)
#-------------------------------------------------------------------------------
# **TODO
# Calculate the force components Fx, Fy, Fz.
mie.gbsc<-function(m,x){
   u<-ofc.abcn(m,x)
   pM<-nrow(u)
   An<-Bn<-Cn<-GTE<-GTM<-l<-m<-j<-1:(pM*(pM+2))
   for(k in 1:pM){
      s<--k:k
      n<-k*(k+1)+s
      l[n]<-k
      m[n]<-s
      j[n]<-n
      GTE[n]<-0
      GTM[n]<-0
      GTE[k*(k+1)+1]<-(1i^k)*sqrt(4*pi*(2*k+1))
      GTM[k*(k+1)+1]<--1i*(1i^k)*sqrt(4*pi*(2*k+1))
      An[n]<-u$An[k]
      Bn[n]<-u$Bn[k]
      Cn[n]<-u$Cn[k]
   }
   return(data.frame(j,l,m,An,Bn,Cn,GTE,GTM))
}
#-------------------------------------------------------------------------------
# SCATTERING COEFFICIENTS
#-------------------------------------------------------------------------------
# Single value at time
mie.qnmx<-function(m,x){
   ab<-mie.anbn(m,x)
   n<-nrow(ab)
   n<-1:n
   an<-ab$an
   bn<-ab$bn
   Qsca<-(2/(x^2))*sum((2*n+1)*(abs(an)^2+abs(bn)^2))
   Qext<-(2/(x^2))*sum((2*n+1)*Re(an+bn))
   Q<-data.frame(Qsca,Qext)
   return(Q)
}
# Vector of values
mie.scat<-function(m,x){
   u<-data.frame(Qsca=x,Qext=x)
   for(i in 1:length(x)){
      u[i,]<-mie.qnmx(m,x[i])
   }
   return(u)
}
#-------------------------------------------------------------------------------
# BESSEL FUNCTIONS AND DERIVATIVES CALCULATED BY DOWNWARD RECURRENCE
#-------------------------------------------------------------------------------
# CYLINDRICAL BESSEL FUNCTIONS [DONE]
cyl.rjyn<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcf.cbld(nmax,x) # Last element
	gn[nmax+1]<-lcf.cbrd(nmax,x) # Last element
	Sn<-(0:nmax)/x
	nj<-(nmax+1):2        # n+1
	Gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n]+Dn[n]
   	  #Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  Gm[n-1]<-Sn[n]+Dm[n]
   	  Dm[n-1]<-Sn[n-1]-1/Gm[n-1]
   	  # modified (permits one step normalization)
   	  gn[n-1]<-Sn[n]*gn[n]+Dn[n]
   	  Dn[n-1]<-Sn[n-1]*gn[n-1]-gn[n]
   	  # Normalization
   	  #print(c(gn[n-1],Dn[n-1]))
   	  if(abs(gn[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(gn[n-1],Dn[n-1]))
   	  	 Dn<-Dn/gn[n-1] # this must be done first
   	  	 gn<-gn/gn[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(gn[n-1],Dn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gn[1])<abs(gn[2])){
      Jn<-(gn/gn[1])*besselJ(x,0) # create functions for normalizations
   }else{
      Jn<-(gn/gn[2])*besselJ(x,1)   	
   }
   # Its Derivative
   if(abs(Dn[1])>abs(Dn[2])){
   	  dJn<-(Dn/Dn[1])*(-Jn[2])
   }else{
   	  dJn<-(Dn/Dn[2])*.5*(Jn[1]-Jn[3])
   }
   #DIRECT CALCULATION
   Jl<-dJl<-Jn
   Jl[1]<-besselJ(x,0)
   dJl[1]<-Jl[1]*Dm[1]
   for(n in 1:(nmax)){
   	  Jl[n+1]<-Jl[n]/Gm[n]
   	  dJl[n+1]<-Jl[n+1]*Dm[n+1]
   }
   # Return results
   return(data.frame(gm=Gm,Dn=Dm,Jn,Jl,dJn,dJl))
}
# SPHERICAL BESSEL FUNCTIONS [DONE]
sph.rjyn<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcf.sbld(nmax,x) # Last element
	gn[nmax+1]<-lcf.sbrd(nmax,x) # Last element
	Sn<-(0:(nmax+1))/x
	nj<-(nmax+1):2        # n+1
	Gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n+1]+Dn[n]
   	  #Dn[n-1]<-Sn[n-1]-1/gn[n-1]
   	  Gm[n-1]<-Sn[n+1]+Dm[n]
   	  Dm[n-1]<-Sn[n-1]-1/Gm[n-1]
   	  # modified (permits one step normalization)
   	  gn[n-1]<-Sn[n+1]*gn[n]+Dn[n]
   	  Dn[n-1]<-Sn[n-1]*gn[n-1]-gn[n]
   	  # Normalization
   	  #print(c(gn[n-1],Dn[n-1]))
   	  if(abs(gn[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(gn[n-1],Dn[n-1]))
   	  	 Dn<-Dn/gn[n-1] # this must be done first
   	  	 gn<-gn/gn[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(gn[n-1],Dn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gn[1])<abs(gn[2])){
      jn<-(gn/gn[1])*sin(x)/x # create functions for normalizations
   }else{
      jn<-(gn/gn[2])*(-cos(x)+sin(x)/x)*(1/x)   	
   }
   # Its Derivative
   if(abs(Dn[1])>abs(Dn[2])){
   	  djn<-(Dn/Dn[1])*(-jn[2])
   }else{
   	  djn<-(Dn/Dn[2])*(1/3)*(jn[1]-2*jn[3])
   }
   # Return results
   return(data.frame(rh=Gm,An=Dm,jn,djn))
}
# RICCATI BESSEL FUNCTIONS [DONE]
ric.rjyn<-function(x,nmax){
	Dn<-rep(0,nmax+1)  # Vector 
	gn<-rep(1,nmax+1)  # Vector
	Dn[nmax+1]<-lcf.rbld(nmax,x) # Last element
	gn[nmax+1]<-lcf.sbrd(nmax,x) # Last element
	Sn<-(0:(nmax))/x
	nj<-(nmax+1):2        # n+1
	Gm<-gn
	Dm<-Dn
	# DOWNWARD RECURRENCE
	RN<-1
   for(n in nj){
   	  # original
   	  #gn[n-1]<-Sn[n]+Dn[n]
   	  #Dn[n-1]<-Sn[n]-1/gn[n-1]
   	  Gm[n-1]<-Sn[n]+Dm[n]
   	  Dm[n-1]<-Sn[n]-1/Gm[n-1]
   	  # modified (permits one step normalization)
   	  gn[n-1]<-Sn[n]*gn[n]+Dn[n]
   	  Dn[n-1]<-Sn[n]*gn[n-1]-gn[n]
   	  # Normalization
   	  #print(c(gn[n-1],Dn[n-1]))
   	  if(abs(gn[n-1])>1e100){
   	  	 cat("renorming...\n")
   	  	 #print(c(gn[n-1],Dn[n-1]))
   	  	 Dn<-Dn/gn[n-1] # this must be done first
   	  	 gn<-gn/gn[n-1] # otherwise the result will be wrong.
   	  }
   	  #print(c(gn[n-1],Dn[n-1]))
   }
   # one step normalization taking care about zeros
   # Bessel function
   if(abs(gn[1])<abs(gn[2])){
      Jn<-(gn/gn[1])*sin(x) # create functions for normalizations
   }else{
      Jn<-(gn/gn[2])*(sin(x)/x-cos(x)) 	
   }
   # Its Derivative
   if(abs(Dn[1])>abs(Dn[2])){
   	  dJn<-(Dn/Dn[1])*(-Jn[2])
   }else{
   	  dJn<-(Dn/Dn[2])*(1/3)*(2*Jn[1]-Jn[3])
   }
   # Return results
   return(data.frame(gm=Gm,Dn=Dm,Jn,dJn))
}
#-------------------------------------------------------------------------------
# TESTING FUNCTIONS
#-------------------------------------------------------------------------------
# Cylindrical Bessel Derivatives for J
tst.cbdj<-function(x,n){
   return(.5*(besselJ(x,n-1)-besselJ(x,n+1)))
}
# Cylindrical Hankel
tst.cbh1<-function(x,n,type=1){
   if(type==1){
      return(besselJ(x,n)+1i*besselY(x,n))
   }
   if(type==2){
      return(besselJ(x,n)-1i*besselY(x,n)) 
   }
}
# Derivative of Cylindrical Hankel
tst.cbdh<-function(x,n,type=1){
   return(.5*(tst.cbh1(x,n-1,type)-tst.cbh1(x,n+1,type)))
}
# Spherical Bessel Functions
tst.sbjn<-function(x,n){
	return(sqrt(pi/2)*besselJ(x,n+.5)/sqrt(x))
}
# Riccati Bessel Functions
tst.rbjn<-function(x,n){
	return(sqrt(pi/2)*besselJ(x,n+.5)*sqrt(x))
}
# Derivative of Spherical Bessel Functions
tst.sbdj<-function(x,n){
   a<-n/(2*n+1)
   b<-(n+1)/(2*n+1)
	return(a*tst.sbjn(x,n-1)-b*tst.sbjn(x,n+1))
}
# Derivative of Riccati Bessel Functions
tst.rbdj<-function(x,n){
   a<-n/(2*n+1)
   b<-(n+1)/(2*n+1)
	return(b*tst.sbjn(x,n-1)-a*tst.sbjn(x,n+1))
}
# Coefficients calculated by definitions
tst.abcd<-function(n,x){
   An<-.5/x+tst.cbdj(x,n+.5)/besselJ(x,n+.5)
   Bn<-.5/x+tst.cbh1(x,n+.5,2)/tst.cbh1(x,n+.5,2)
   Cn<-besselJ(x,n+.5)/tst.cbh1(x,n+.5,2)
   return(data.frame(An,Bn,Cn))
}
tst.anbn<-function(n,m,x){
   u.1<-tst.abcd(n,x)
   u.m<-tst.abcd(n,m*x)
   Ta<-(u.m$An/m-u.1$An)/(u.m$An/m-u.1$Bn)
   Tb<-(u.m$An*m-u.1$An)/(u.m$An*m-u.1$Bn)
   an<-u.1$Cn*Ta
   bn<-u.1$Cn*Tb
   return(data.frame(An=u.1$An,Bn=u.1$Bn,Cn=u.1$Cn,Ta,Tb,an,bn))
}





#-------------------------------------------------------------------------------
# COMPARISONS
#-------------------------------------------------------------------------------
# RATIO BETWEEN SPHERICAL OR RICCATI J/Y???
cmp.rhcj<-function(n,x){
   a<-besselJ(x,n+.5)
   b<-besselY(x,n+.5)
   return(1/(1+1i*b/a))
}
# RHO FOR SPHERICAL AND RICCATI BESSEL J
cmp.rhsj<-function(n,x){
   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
   b<-besselJ(x,n+1.5)/besselJ(x,n+2.5)
   c<-(2*n+3)/x
   cat("rho[n]=j[n]/j[n+1]\n")
   cat("a=rho[n]\n")
   cat("b=rho[n+1]\n")
   cat("c=(2n+3)/x\n")
   cat("a+1/b=c\n")
   return(data.frame(a,b,u=a+1/b,c))
}
# RHO FOR SPHERICAL AND RICCATI BESSEL Y
cmp.rhsy<-function(n,x){
   a<-besselY(x,n+ .5)/besselY(x,n+1.5)
   b<-besselY(x,n+1.5)/besselY(x,n+2.5)
   c<-(2*n+3)/x
   cat("rho[n]=y[n]/y[n+1]\n")
   cat("a=rho[n]\n")
   cat("b=rho[n+1]\n")
   cat("c=(2n+3)/x\n")
   cat("a+1/b=c\n")
   return(data.frame(a,b,u=a+1/b,c))
}
# RHO FOR SPHERICAL AND RICCATI BESSEL H
cmp.rhsh<-function(n,x,k=1){
   a<-(besselJ(x,n+ .5)+k*1i*besselY(x,n+ .5))/
      (besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))
   b<-(besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))/
      (besselJ(x,n+2.5)+k*1i*besselY(x,n+2.5))
   c<-(2*n+3)/x
   cat("rho[n]=h[n]/h[n+1]\n")
   cat("a=rho[n]\n")
   cat("b=rho[n+1]\n")
   cat("c=(2n+3)/x\n")
   cat("a+1/b=c\n")
   return(data.frame(a,b,u=a+1/b,c))
}


#-------------------------------------------------------------------------------
# CALCULATION TESTS
#-------------------------------------------------------------------------------
# WAVELENGTH (MICROMETERS)
l<-.532
# INDEX OF REFRACTION POLIESTIRENE
np<-1.5725+.0031080/l^2+.00034779/l^4
# INDEX OF REFRACTION POLIESTIRENE
nw<-1.3253+.0027553/l^2+.00003779/l^4
m<-np/nw