#-------------------------------------------------------------------------------
ZJ <-matrix(c( 2.4048, 3.8317, 5.1356, 6.3802, 7.5883, 8.7715,
               5.5201, 7.0156, 8.4172, 9.7610,11.0647,12.3386,
               8.6537,10.1735,11.6198,13.0152,14.3725,15.7002,
              11.7915,13.3237,14.7960,16.2235,17.6160,18.9801,
              14.9309,16.4706,17.9598,19.4094,20.8269,22.2178),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
ZdJ<-matrix(c( 3.8317, 1.8412, 3.0542, 4.2012, 5.3175, 6.4156,
               7.0156, 5.3314, 6.7061, 8.0152, 9.2824,10.5199,
              10.1735, 8.5363, 9.9695,11.3459,12.6819,13.9872,
              13.3237,11.7060,13.1704,14.5858,15.9641,17.3128,
              16.4706,14.8636,16.3475,17.7887,19.1960,20.5755),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
lmax<-50
#-------------------------------------------------------------------------------
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP=200
R<-5*lambda
M<-5  # 0<=M<=5
N<-3  # 1<=N<=5
TM<-FALSE
S=-1
#
x<-seq(-R,R,by=2*R/(NP-1))
y<-x
z<-x
#-------------------------------------------------------------------------------
# MODOS: M <- M+1 because J_0(x) ocupies column 1
# TE - ZdJ
# TM - ZJ
#-------------------------------------------------------------------------------
if(TM){
   g<-ZJ[N,M+1]/R
   MODE<--1
   MODE.PWE<-4
}else{
   g<-ZdJ[N,M+1]/R
   MODE<-1
   MODE.PWE<-3
}
kz<-sqrt(k^2-g^2)
xo<-sample(x,1)
yo<-sample(y,1)
zo<-sample(z,1)
xo<-yo<-zo<-0
x<-x-xo
y<-y-yo
z<-z-zo
xe<-sample(x,1)
ye<-sample(y,1)
ze<-sample(z,1)
if(sqrt(xe^2+ye^2+ze^2)>R){
   cat("::::::::::::::::::::::::::::: r>R :::::::::::::::::::::::::::::\n")}
#
E.TM.H.TE.m<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M-S,s=S)*( 1i*S*kz/g)/sqrt(2)
E.TM.H.TE.z<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M,s=S)
E.TM.H.TE.p<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M+S,s=S)*(-1i*S*kz/g)/sqrt(2)
#
H.TM.E.TE.m<--vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M-S,s=S)*(k/g)*S/sqrt(2)
H.TM.E.TE.z<-0                                              
H.TM.E.TE.p<--vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M+S,s=S)*(k/g)*S/sqrt(2)
#
E.TM.H.TE<-c(E.TM.H.TE.m,E.TM.H.TE.z,E.TM.H.TE.p)
H.TM.E.TE<-MODE*c(H.TM.E.TE.m,H.TM.E.TE.z,H.TM.E.TE.p)
#-------------------------------------------------------------------------------
if(TM){
   E.wfd<-c(E.TM.H.TE.m,E.TM.H.TE.z,E.TM.H.TE.p)
   H.wfd<-c(H.TM.E.TE.m,H.TM.E.TE.z,H.TM.E.TE.p)
}else{
   E.wfd<-c(H.TM.E.TE.m,H.TM.E.TE.z,H.TM.E.TE.p)
   H.wfd<-c(E.TM.H.TE.m,E.TM.H.TE.z,E.TM.H.TE.p)
}
#-------------------------------------------------------------------------------
t<-vwfd.cwg(!TM,M,S,g,kz,xe+xo,ye+yo,ze+zo)
u<-vswf.cwg(TM,M,S,g,kz,xo,yo,zo,lmax)
v<-vswf.hmp(k,xe,ye,ze,lmax)
w<-vswf.pwe(k,xe,ye,ze,lmax,u$GTE,u$GTM)
#-------------------------------------------------------------------------------
Em.hmp<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.hmp<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.hmp<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.hmp<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.hmp<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.hmp<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#-------------------------------------------------------------------------------
E.vwf<-c(t$Em,t$Ez,t$Ep)
H.vwf<-c(t$Hm,t$Hz,t$Hp)

E.hmp<-c(Em.hmp,Ez.hmp,Ep.hmp)
H.hmp<-c(Hm.hmp,Hz.hmp,Hp.hmp)

E.pwe<-c(w$Em,w$Ez,w$Ep)
H.pwe<-c(w$Hm,w$Hz,w$Hp)
#-------------------------------------------------------------------------------
VWF<-as.data.frame(cbind(E.vwf,H.vwf),row.names=c("m","z","p"))
WFD<-as.data.frame(cbind(E.wfd,H.wfd),row.names=c("m","z","p"))
HMP<-as.data.frame(cbind(E.hmp,H.hmp),row.names=c("m","z","p"))
PWE<-as.data.frame(cbind(E.pwe,H.pwe),row.names=c("m","z","p"))
#
print(VWF)
print(WFD)
print(HMP)
print(PWE)
#-------------------------------------------------------------------------------