#-------------------------------------------------------------------------------
# GENERIC WAVE GUIDES
#-------------------------------------------------------------------------------
vswf.gwg<-function(gama,kz,lmax){
   k<-sqrt(gama^2+kz^2)
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-vswf.qlm(kz/k,lmax)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   ll<-u$l
   mm<-u$m
   llp1<-1/sqrt(ll*(ll+1))
   llp1[1]<-0
   #----------------------------------------
   A<-2*(1i^ll)*((k/gama)^2)*Qlm*mm*llp1
   B<-2*(1i^(ll-1))*dQlm*llp1
   return(data.frame(A,B))
}
