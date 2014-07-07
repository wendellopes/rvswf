#-------------------------------------------------------------------------------
# SCATTERING COEFFICIENTS
#-------------------------------------------------------------------------------
# Single value at time
lmie.cal<-function(m,x){
   ab<-lmie.exp(m,x)
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
lmie.sct<-function(m,x){
   u<-data.frame(Qsca=x,Qext=x)
   for(i in 1:length(x)){
      u[i,]<-lmie.cal(m,x[i])
   }
   return(u)
}
