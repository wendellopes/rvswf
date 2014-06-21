#-------------------------------------------------------------------------------
# RECTANGULAR WAVE GUIDE
#-------------------------------------------------------------------------------
vwfd.rwg<-function(TE=TRUE,kx,ky,kz,x,y,z){
if(TE)
{
   te<-1
}else{
   te<-0
}
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("vwfd_rwg",
   TE=as.integer(te),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   kx=as.double(kx) ,ky=as.double(ky) ,kz=as.double(kz),
   x=as.double(x) ,y=as.double(y) ,z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
