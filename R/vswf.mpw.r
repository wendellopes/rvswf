#-------------------------------------------------------------------------------
# MIE PLANE WAVE (PROPAGATING IN Z DIRECTION)
#-------------------------------------------------------------------------------
vswf.mpw<-function(lmax,norm=FALSE,s=1){
#-------------------------------------------------------------------------------
   if(lmax<1){lmax<-1}                         # Pelo menos 1 termo
   LMAX=lmax*(lmax+2)+1                        # Vetor para lmax
#-------------------------------------------------------------------------------
   gte<-rep(0,LMAX)
   gtm<-rep(0,LMAX)
   glm<-function(l,m,norm=FALSE){
      k<-(1i^l)*sqrt(4*pi*(2*l+1))
      k<-rep(k,2*l+1)
      k[m!=s]<-0
      if(norm){
         k<-k/sqrt(2)
      }
      return(k)
   }
   for(l in 0:lmax){
      m<--l:l
      gte[vswf.jlm(l,m)]<-glm(l,m,norm)
   }
   gtm<--1i*s*gte
   U<-data.frame(GTE=gte,GTM=gtm)
   return(U)
}
