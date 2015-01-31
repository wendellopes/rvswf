# Philip Laven Mie Plot Page
lmax<-50
m<-1.33
x<-2*pi*10/.65
# Wriedt Book, pg 65
x<-2*pi*.7/.5
m<-1.5
m<-1.5+.1i
#
th<-seq(0,pi,pi/1999)
u<-lmie.phs(th,m,x)
i.p.o<-abs(u$S2)**2 # paralell
i.s.o<-abs(u$S1)**2 # senkrecht (orthogonal)
i.t.o<-.5*(i.p.o+i.s.o)
i.p<-log10(i.p.o)
i.s<-log10(i.s.o)
i.t<-log10(i.t.o)
#P<-(i.s-i.p)/(i.s+i.p)
# Range
imo<-ceiling(max(abs(c(i.p.o,i.s.o,i.t.o))))
imx<-ceiling(max(abs(c(i.p,i.s,i.t))))-1
# POLAR PLOT
plot(cos(th)*i.s,sin(th)*i.s,type='l',lwd=1,xlab="",ylab="",asp=1,ylim=imx*c(-1,1))
points(cos(-th)*i.p,sin(-th)*i.p,type='l',lwd=1,col='red')
points(cos(th)*i.t,sin(th)*i.t,type='l',lwd=1,col='blue')
points(cos(-th)*i.t,sin(-th)*i.t,type='l',lwd=1,col='blue')
ph<-seq(0,2*pi,2*pi/199)
for(i in 1:imx){
   points(i*cos(ph),i*sin(ph),type='l',col='lightblue')
}
for(i in 1:36){
   points(c(0,imx*cos(pi*i/18)),c(0,imx*sin(pi*i/18)),type='l',col='lightblue')
}
# LINEAR PLOT
plot(th*180/pi,i.s.o/imo,type='l',log='y',ylim=c(1e-5,1e0))
points(th*180/pi,i.p.o/imo,col='red',type='l')
points(th*180/pi,i.t.o/imo,col='blue',type='l')
