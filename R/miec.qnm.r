#-------------------------------------------------------------------------------
# SCATTERING COEFFICIENTS
#-------------------------------------------------------------------------------
# Single value at time
mie.qnmx<-function(m,x){
   ab<-mie.anbn(m,x)
   n<-nrow(ab)
   n<-1:n
   an<-ab$an
   bn<-ab$bn
   Qsca<-(2/(x^2))*sum((2*n+1)*(abs(an)^2+abs(bn)^2))
   Qext<-(2/(x^2))*sum((2*n+1)*Re(an+bn))
   Q<-data.frame(Qsca,Qext)
   return(Q)
}
# Vector of values
mie.scat<-function(m,x){
   u<-data.frame(Qsca=x,Qext=x)
   for(i in 1:length(x)){
      u[i,]<-mie.qnmx(m,x[i])
   }
   return(u)
}
