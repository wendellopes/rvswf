#-------------------------------------------------------------------------------
# J_{n+1}/J_{n} [OK] BARNETT - Cylindrical Bessel Ratio Inverse
#-------------------------------------------------------------------------------
lcfe.cbi<-function(n,x,NMAX=2000,code="C"){
   nmaxo<-NMAX
   fn<-0
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      u<-.C("lcfe_cbi",
         n=as.integer(n),
         x=as.double(x),
         NMAX=as.integer(NMAX),
         fn=as.double(fn)
      )
      if(u$NMAX>nmaxo){
         stop(paste("The accuracy was not reached in",NMAX,"iterations!"))
      }
      return(u$fn)
   }
   else{
      n<-as.integer(n)
      # Constants
      eo<-.Machine$double.xmin
      ACC<-10^-50
      # initialization of calculations
      fn<-0    # bo
      if(abs(fn)<eo){fn<-eo} # migth be zero
      Pn<-fn
      Qn<-0
      # Loop Parameters
      j<-0
      Dn<-10
      while(abs(Dn-1)>ACC){
         if(j>NMAX){
            NMAX<-j
            break
         }
         j<-j+1
         an<-(-1)^sign(j-1) # j=1, an=1; j>1, an=-1
         bn<-lcfe.afs(2*(n+j),x);
         Pn<-bn+an/Pn
         if(abs(Pn)<eo){Pn<-eo} # migth be zero
         Qn<-bn+an*Qn
         if(abs(Qn)<eo){Qn<-eo} # migth be zero
         Qn<-1/Qn
         Dn<-Pn*Qn
         fn<-fn*Dn
      }
      if(NMAX>nmaxo){
         stop(paste("The accuracy was not reached in",NMAX,"iterations!"))
      }
   return(fn)
   }
}
