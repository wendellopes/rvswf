#-------------------------------------------------------------------------------
# NAMES
#-------------------------------------------------------------------------------
md<-ifelse(TM,"tm","te")
nm.pwe<-ifelse(lmax<10, paste("bbz.vswf.",md,".0",lmax,sep=""),
                        paste("bbz.vswf.",md,"." ,lmax,sep=""))
nm.vwf<-paste("bbz.vwfd.",md,".00",sep="")
nm.vwf.i<-paste(nm.vwf,".png",sep="")
nm.vwf.d<-paste(nm.vwf,".Rdata",sep="")
nm.pwe.i<-paste(nm.pwe,".png",sep="")
nm.pwe.d<-paste(nm.pwe,".Rdata",sep="")
#-------------------------------------------------------------------------------
# IMAGE
#-------------------------------------------------------------------------------
source("Image.R")
u<-seq(0,2*pi,pi/200)
if(TM){
   zl<-range(Re(tez.BBZ))
   #1
   if(!file.exists(nm.vwf.i)){
      png(nm.vwf.i)
      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.BBZ),nlevels=256,axes=TRUE,
         color.palette=cm.colors,asp=1,#zlim=zl,
         plot.axes={
            axis(1);
            axis(2);
            abline(h=xo/lambda,v=yo/lambda,col='green')
         },
         xlab=expression(y/lambda),ylab=expression(x/lambda)
      )
      dev.off()
   }
   #2
   png(nm.pwe.i)
      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.PWE),nlevels=256,axes=TRUE,
         color.palette=cm.colors,asp=1,#zlim=zl,
         plot.axes={
            axis(1);
            axis(2);
            abline(h=xo/lambda,v=yo/lambda,col='green')
         },
         xlab=expression(y/lambda),ylab=expression(x/lambda)
   )
   dev.off()
}else{
   zl<-range(Re(thz.BBZ))
   #1
   if(!file.exists(nm.vwf.i)){
      png(nm.vwf.i)
      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.BBZ),nlevels=256,axes=TRUE,
         color.palette=cm.colors,asp=1,#zlim=zl,
         plot.axes={
            axis(1);
            axis(2);
            abline(h=xo/lambda,v=yo/lambda,col='green')
            },
         xlab=expression(y/lambda),ylab=expression(x/lambda)
      )
      dev.off()
   }
   #2
   png(nm.pwe.i)
      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.PWE),nlevels=256,axes=TRUE,
         color.palette=cm.colors,asp=1,#zlim=zl,
         plot.axes={
            axis(1);
            axis(2);
            abline(h=xo/lambda,v=yo/lambda,col='green')
         },
         xlab=expression(y/lambda),ylab=expression(x/lambda)
   )
   dev.off()
}
#-------------------------------------------------------------------------------
# DATASETS
#-------------------------------------------------------------------------------
if(TM){
   save(BBZ, file=nm.vwf.d)
   save(PWE, file=nm.pwe.d)
}else{
   save(BBZ, file=nm.vwf.d)
   save(PWE, file=nm.pwe.d)
}
