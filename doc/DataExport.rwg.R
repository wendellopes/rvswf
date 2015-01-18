#-------------------------------------------------------------------------------
# SINGLE LINE PLOT - ALONG y DIRECTION AT x=xo
#-------------------------------------------------------------------------------
if(!TM){
   RWG<-vwfd.rwg(TE=!TM,kx,ky,kz,x+xo,yo,zo)$Hz
}else{
   RWG<-vwfd.rwg(TE=!TM,kx,ky,kz,x+xo,yo,zo)$Ez
}
ii<-seq(2,32,1)
ss<-c(2,5,11,17,23,31)
#ss<-c(5,11,23,31)
cores<-hsv(h=ii/40)
VPWE<-array(0,c(NPX,1+length(ii)))
VPWE[,1]<-RWG
tlg<-array(0,length(ii))
j<-1
for(i in ii){
   j<-j+1
   u<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,i)
   if(!TM){
      VPWE[,j]<-vswf.pwe(k,x,0,0,i,u$GTE,u$GTM)$Hz
   }else{
      VPWE[,j]<-vswf.pwe(k,x,0,0,i,u$GTE,u$GTM)$Ez
   }
}
#-------------------------------------------------------------------------------
VPWE[is.na(VPWE)]<-0 
#-------------------------------------------------------------------------------
# Exporting File
RPWE<-as.data.frame(apply(VPWE,2,Re))
IPWE<-as.data.frame(apply(VPWE,2,Im))
names(RPWE)<-c("EX",paste("L",ii,sep=""))
names(IPWE)<-c("EX",paste("L",ii,sep=""))
RPWE$X<-IPWE$X<-(x+xo)/lambda
write.table(file="rpwe.dat",x=RPWE,row.names=FALSE,col.names=TRUE)
write.table(file="ipwe.dat",x=IPWE,row.names=FALSE,col.names=TRUE)
#-------------------------------------------------------------------------------
# saving the files
#-------------------------------------------------------------------------------
if(TM){
   save(RWG, file="RWG.TM.Rdata")
   save(PWE, file="PWE.TM.Rdata")
}else{
   save(RWG, file="RWG.TE.Rdata")
   save(PWE, file="PWE.TE.Rdata")
}
