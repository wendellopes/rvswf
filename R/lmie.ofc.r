#' Lorentz-Mie Optical Force Calculation.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param by \code{LD} using logarithmic derivative; \code{RB} using the ratio
#' between Ricatti-Bessel functions.
#' @return The coefficients \eqn{A_n}, \eqn{B_n} and \eqn{C_n} used to calculate
#' the Optical force.
#' @include lmie.exp
#' @export
#' @examples
#' lmax<-20
#' lambda=.532  # Green 532 nm
#' m<-ridx.pol(lambda)/ridx.h2o(lambda)
#' g<-vswf.mpw(lmax+1)
#' g<-vswf.bbp(3,4,1,2,3,lmax+1,1,1,1)
#' x<-20
#' 
#' F<-lmie.ofc(m,x,g$GTE,g$GTM,lmax)
#' print(F)
lmie.ofc<-function(m,x,GTE,GTM,lmax=floor(x+2+4*x^(1/3)),code="C"){
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,lmax)
      Fx<-Fy<-Fz<-0
      u<-.C("lmie_ofc",
         M=as.complex(m),           # m
         x=as.complex(x),           # x
         lmax=as.integer(lmax),     # lmax
         NMAX=as.integer(200),      # NMAX
         GTE=as.complex(GTE),
         GTM=as.complex(GTM),
         Fx=as.double(Fx),
         Fy=as.double(Fy),
         Fz=as.double(Fz))    
      return(data.frame(Fx=u$Fx,Fy=u$Fy,Fz=u$Fz))
   }else{
      u<-lmie.exp(m,x,by="LD",lmax=lmax+1)   
      #print(u)
      nu<-lmax
      nm<-nu*(nu+2)
      Anm<-Bnm<-Cnm<-1:nm
      # Calculations
      n0<-1:nu
      n1<-n0+1
      An<-u$an[n0]+Conj(u$an[n1])-2*u$an[n0]*Conj(u$an[n1])
      Bn<-u$bn[n0]+Conj(u$bn[n1])-2*u$bn[n0]*Conj(u$bn[n1])
      Cn<-u$an[n0]+Conj(u$bn[n0])-2*u$an[n0]*Conj(u$bn[n0])
      #print(data.frame(An,Bn,Cn))
      #print(data.frame(GTE,GTM))
      # Results
      Fx<-0
      Fy<-0
      Fz<-0
      for(l in 1:lmax){
         for(m in -l:l){
            K1<-sqrt((l*(l+2))/((2*l+1)*(2*l+3)))
            K2<-sqrt((l+m+2)*(l+m+1))
            K3<-sqrt((l-m)*(l+m+1))
            K4<-sqrt((l+m+1)*(l-m+1))
   
            U<-An[l]*GTM[vswf.jlm(l,m)]*Conj(GTM[vswf.jlm(l+1,m+1)])+
               Bn[l]*GTE[vswf.jlm(l,m)]*Conj(GTE[vswf.jlm(l+1,m+1)])+
               Conj(An[l]*GTM[vswf.jlm(l,-m)])*GTM[vswf.jlm(l+1,-m-1)]+
               Conj(Bn[l]*GTE[vswf.jlm(l,-m)])*GTE[vswf.jlm(l+1,-m-1)]
            V<-Cn[l]*GTM[vswf.jlm(l,m)]*Conj(GTE[vswf.jlm(l,m+1)])-
                 Conj(Cn[l])*GTE[vswf.jlm(l,m)]*Conj(GTM[vswf.jlm(l,m+1)])
            W<-An[l]*GTM[vswf.jlm(l,m)]*Conj(GTM[vswf.jlm(l+1,m)])-
                 Bn[l]*GTE[vswf.jlm(l,m)]*Conj(GTE[vswf.jlm(l+1,m)])
            Y<-Cn[l]*GTM[vswf.jlm(l,m)]*Conj(GTE[vswf.jlm(l,m)])
   
            Fxy<-(K1*K2*U-K3*V/l)/(l+1)
            Fzz<-(K1*K4*W-m*Y/l)/(l+1)
         
            #print(c(Fxy,Fzz))

            Fx<-Fx+Re(Fxy)
            Fy<-Fy+Im(Fxy)
            Fz<-Fz+Re(1i*Fzz)
         }
      }
      return(data.frame(Fx,Fy,Fz))
   }
}
