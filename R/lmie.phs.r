#' Phase calculation
#' 
#' @param th The angle at which the phase will be calculated
#' @param m  The relative refraction index
#' @param x  The form factor
#' @param lmax The maximum value of \code{n}
#' @examples
#' # Bohren and Hoffman example, pg 114, section 4.4.5
#' th<-seq(0,pi,pi/199)
#' lmax<-5
#' m<-1.33+1e-8i
#' x<-3
#' u<-lmie.phs(th,m,x,lmax)
#' i.p<-abs(u$S2)**2 # paralell
#' i.s<-abs(u$S1)**2 # senkrecht (orthogonal)
#' P<-(i.s-i.p)/(i.s+i.p)
#' # Range
#' imx<-max(c(i.p,i.s))
#' imn<-min(c(i.p,i.s))
#' # Main plot
#' plot(cos(th)*i.s,sin(th)*i.s,type='l',lwd=1,asp=1,xlim=imx*c(-1,1),xlab="",ylab="")
#' points(cos(-th)*i.p,sin(-th)*i.p,type='l',lwd=1,col='red')
#' # x10 plot (zoom)
#' points(10*cos(th)*i.s,10*sin(th)*i.s,type='l',lwd=2)
#' points(10*cos(-th)*i.p,10*sin(-th)*i.p,type='l',lwd=2,col='red')
#' # Log plot
#' plot(th,i.p,type='l',log='y',ylim=c(imn,imx))
#' points(th,i.s,type='l',col='red')
#' # Polarization
#' plot(th,P,type='l')
lmie.phs<-function(th,m,x,lmax=floor(x+2+3*x^{1/3})){
   if(length(th)>1){
      SS<-sapply(th,lmie.phs,m=m,x=x,lmax=lmax)
      SS<-as.data.frame(apply(t(SS),2,as.complex))
      return(SS)
   }else{
      u<-lmie.exp(m,x,NMAX=lmax)
      p<-t<-rep(0,lmax)
      w<-cos(th)
      # Initial value
      #p[0]<-0           # \pi_0
      p[1]<-1            # \pi_1
      p[2]<-3*w          # \pi_2
      #t[0]<-0           # \tau_0
      t[1]<-w            # \tau_1
      t[2]<-3*(2*w**2-1) # 3*cos(2*th)  # \tau_2
      for(j in 3:lmax){# j=n+1
         p[j]<-((2*j-1)/(j-1))*w*p[j-1]-(j/(j-1))*p[j-2]
         t[j]<-j*w*p[j]-(j+1)*p[j-1]
      }
      n<-1:lmax
      S1<-sum(((2*n+1)/(n*(n+1)))*(u$an*p+u$bn*t))
      S2<-sum(((2*n+1)/(n*(n+1)))*(u$an*t+u$bn*p))
      return(data.frame(S1,S2))
   }
}