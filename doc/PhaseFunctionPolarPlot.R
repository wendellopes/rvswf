th<-seq(0,pi,pi/1999)
lmax<-50
m<-1.33
x<-2*pi*10/.65
u<-lmie.phs(th,m,x)
i.p<-abs(u$S2)**2 # paralell
i.s<-abs(u$S1)**2 # senkrecht (orthogonal)
i.t<-log10(.5*(i.p+i.s))
i.p<-log10(i.p)
i.s<-log10(i.s)
#P<-(i.s-i.p)/(i.s+i.p)
# Range
imx<-ceiling(max(abs(c(i.p,i.s,i.t))))-1
# Main plot
#plot(cos(th)*i.t,sin(th)*i.t,type='l',lwd=1,xlab="",ylab="")
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
