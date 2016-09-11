#' 
#' n<-1.1;thmax<-79
#' n<-1.2;thmax<-72
#' #n<-1.4;thmax<-64
#' #n<-1.6;thmax<-60
#' #n<-1.8;thmax<-59
#' #n<-2.0;thmax<-59
#' #n<-2.5;thmax<-64
#' 
#' th<-seq(0,90,.1)*pi/180
#' a<-lmie.rfc(th,1/n,pp=0)
#' b<-lmie.rfc(th,1/n,pp=1)
#' 
#' tg<-th*180/pi
#' plot(tg,a$Ft,type='l',xlab=expression(theta),ylab="Q",ylim=range(c(a$Ft,b$Ft)))
#' points(tg,abs(a$Fg),type='l',col='blue')
#' points(tg,abs(a$Fs),type='l',col='red')
#' 
#' points(tg,abs(b$Ft),lty=2,type='l')
#' points(tg,abs(b$Fg),lty=2,type='l',col='blue')
#' points(tg,abs(b$Fs),lty=2,type='l',col='red')
#' 
#' legend("topleft",c("Full S","Grad S","Scat S","Full P","Grad P","Scat P"),
#'       lty=c(1,1,1,2,2,2),col=c("black","blue","red"))
#' 
#' #abline(v=thmax,col='green')
#' #print(round(min(a$Fg),3))
#' 
lmie.rfc<-function(th.i,n,pp=.5){
   # Snell law
   #
   # n1*sin th1 = n2*sin th2 -> sin th2 = n sin th1
   # n = n1/n2
   #
   # Fresnel coefficients
   #
   # Rs=((n1*cos thi-n2*cos thr)/(n1*cos thi + n2*cos thr))**2
   # Rp=((n1*cos thr-n2*cos thi)/(n1*cos thr + n2*cos thi))**2
   #
   # Rs=((n*cos thi-cos thr)/(n*cos thi + cos thr))**2
   # Rp=((n*cos thr-cos thi)/(n*cos thr + cos thi))**2
   #
   th.r<-asin(sin(th.i)*n)
   RS<-((n*cos(th.i)-cos(th.r))/(n*cos(th.i)+cos(th.r)))**2
   RP<-((n*cos(th.r)-cos(th.i))/(n*cos(th.r)+cos(th.i)))**2
   R<-pp*RP+(1-pp)*RS
   T<-1-R
   Fo<-T**2/(1+R**2+2*R*cos(2*th.r))
   Fs<-1+R*cos(2*th.i)-(cos(2*th.i-2*th.r)+R*cos(2*th.i))*Fo
   Fg<-  R*sin(2*th.i)-(sin(2*th.i-2*th.r)+R*sin(2*th.i))*Fo
   Ft<-sqrt(Fs**2+Fg**2)
   return(data.frame(Ft,Fs,Fg))
} 
