#-------------------------------------------------------------------------------
# PARAMETERS
#-------------------------------------------------------------------------------
TM<-TRUE
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP=200
a<-5*lambda
b<-8*lambda
M<-7
N<-6
kx<-M*pi/a
ky<-N*pi/b
gamma<-sqrt(kx^2+ky^2)
kz<-sqrt(k^2-gamma^2)
dx<-a/(NP-1)
dy<-b/(NP-1)
x<-seq(0,a,dx)
y<-seq(0,b,dy)
z<-x
lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
#-------------------------------------------------------------------------------
nxo<-sample(1:200,1)
nyo<-sample(1:200,1)
nzo<-sample(1:200,1)
# insere desvio em (xo,yo,zo) para evitar NaN
xo<-x[nxo]+dx/2
yo<-y[nyo]+dy/2
zo<-z[nzo]+dx/2
#
xo<-x[27]+dx/2
yo<-y[82]+dy/2
zo<-z[29]+dx/2
#xo<-0
#yo<-0
#zo<-0
#-------------------------------------------------------------------------------
# MUDA O REFERENCIAL
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
# POSICAO DO CAMPO TESTADO (DE 1 A 200)
#-------------------------------------------------------------------------------
nx<-sample(1:200,1)
ny<-sample(1:200,1)
nz<-sample(1:200,1)
xe<-x[nx]
ye<-y[ny]
ze<-z[nz]
xe<-x[73]
ye<-y[80]
ze<-z[173]
#xe<-5e-7
#ye<-0
#ze<-xe
#-------------------------------------------------------------------------------
# Calcula os campos
#-------------------------------------------------------------------------------
# MODO TRANSVERSAL MAGNETICO
E.TM.x<-1i*((kz*kx)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.y<-1i*((kz*ky)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.z<-sin(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
H.TM.x<--1i*((k*ky)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.y<- 1i*((k*kx)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.z<- 0
#
E.TM.m<-(E.TM.x-1i*E.TM.y)/sqrt(2)
E.TM.p<-(E.TM.x+1i*E.TM.y)/sqrt(2)
E.TM<-c(E.TM.m,E.TM.z,E.TM.p)
#
H.TM.m<-(H.TM.x-1i*H.TM.y)/sqrt(2)
H.TM.p<-(H.TM.x+1i*H.TM.y)/sqrt(2)
H.TM<-c(H.TM.m,H.TM.z,H.TM.p)
#-------------------------------------------------------------------------------
# MODO TRANSVERSAL ELETRICO
#-------------------------------------------------------------------------------
H.TE.x<--1i*((kz*kx)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.y<--1i*((kz*ky)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.z<- cos(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
E.TE.x<--1i*((k*ky)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.y<- 1i*((k*kx)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.z<- 0
#
#
E.TE.m<-(E.TE.x-1i*E.TE.y)/sqrt(2)
E.TE.p<-(E.TE.x+1i*E.TE.y)/sqrt(2)
E.TE<-c(E.TE.m,E.TE.z,E.TE.p)
#
H.TE.m<-(H.TE.x-1i*H.TE.y)/sqrt(2)
H.TE.p<-(H.TE.x+1i*H.TE.y)/sqrt(2)
H.TE<-c(H.TE.m,H.TE.z,H.TE.p)
#-------------------------------------------------------------------------------
if(TM){
   E.wfd<-E.TM
   H.wfd<-H.TM
}else{
   E.wfd<-E.TE
   H.wfd<-H.TE
}
#-------------------------------------------------------------------------------
t<-vwfd.rwg(!TM,kx,ky,kz,xe+xo,ye+yo,ze+zo)
u<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,lmax)
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
