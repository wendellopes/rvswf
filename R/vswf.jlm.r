#-------------------------------------------------------------------------------
# POSITIONER
#-------------------------------------------------------------------------------
vswf.jlm<-function(l,m){
   j.lm<-l*(l+1)+m+1    
   j.tf<-abs(m)>l
   j.lm[j.tf]<-1
   return(j.lm)
}