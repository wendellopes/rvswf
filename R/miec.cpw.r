#-------------------------------------------------------------------------------
# G^TE, G^TM For Plane Waves +1 (Jackson's Classical Electrodynamics book)
#-------------------------------------------------------------------------------
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
