#-------------------------------------------------------------------------------
# ESCREVER PARA CONVERTER LOCAL
SEP<-"#-------------------------------------------------------------------------------\n"
#-------------------------------------------------------------------------------
source("BSC.R")
source("HansenMultipoles.R")
#-------------------------------------------------------------------------------
# PARAMETROS DO FEIXE DE BESSEL
#-------------------------------------------------------------------------------
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP<-200
a<-5*lambda
b<-8*lambda
M<-3
N<-5
kx<-M*pi/a
ky<-N*pi/b
gama<-sqrt(kx^2+ky^2)
kz<-sqrt(k^2-gama^2)
x<-seq(0,a,by=a/(NP-1))
y<-seq(0,b,by=b/(NP-1))
z<-x
lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
cat(SEP)
cat("LMAX=",lmax,"\n",sep="")
#-------------------------------------------------------------------------------
xo<-sample(x,1)
yo<-sample(y,1)
zo<-sample(z,1)
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
xe<-sample(x,1)
ye<-sample(y,1)
ze<-sample(z,1)
#-------------------------------------------------------------------------------
S<--1
M<-1
#-------------------------------------------------------------------------------
# MODO TE
# MODO TM -> E <- H, H <- -E
#-------------------------------------------------------------------------------
Em<- 1i*S*((kz/k)*(gama/k))*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Ez<-((gama/k)^2)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M,S)
Ep<--1i*S*((kz/k)*(gama/k))*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#
Hm<-S*(gama/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Hz<-0
Hp<-S*(gama/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#-------------------------------------------------------------------------------
E<-c(Em,Ez,Ep)
H<-c(Hm,Hz,Hp)
cat(SEP)
print(cbind(E,H))
cat(SEP)
#-------------------------------------------------------------------------------
u<-BesselBeamTE(gama,kz,xo,yo,zo,lmax,M,S)
v<-HansenMultipoles(k,xe,ye,ze,lmax)

Em.pwe<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.pwe<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.pwe<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.pwe<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.pwe<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.pwe<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#-------------------------------------------------------------------------------
# CONFERINDO
#-------------------------------------------------------------------------------
E.pwe<-c(Em.pwe,Ez.pwe,Ep.pwe)
H.pwe<-c(Hm.pwe,Hz.pwe,Hp.pwe)
cat(SEP)
print(cbind(E.pwe,H.pwe))
cat(SEP)
#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
LE<-FALSE
if(LE){
   source("LE-BBZ.R")
}
#-------------------------------------------------------------------------------
# FIM
#-------------------------------------------------------------------------------
