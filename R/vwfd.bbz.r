#-------------------------------------------------------------------------------
# BESSEL BEAM Z
#-------------------------------------------------------------------------------
vwfd.bbz<-function(TM=TRUE,M,S,g,kz,x,y,z){
   if(abs(S)!=1){
      stop("S must be plus or minus 1")
   }
if(TM){
   tm<-0
}else{
   tm<-1
}
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("vwfd_bbz",
   TM=as.integer(tm),M=as.integer(M),S=as.integer(S),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   gama=as.double(g),kz=as.double(kz),
   x=as.double(x),y=as.double(y),z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
