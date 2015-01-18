library(gdata)
#for(k in c("cwg","rwg","bbp","bbz")){
for(k in c("bbp")){
   for(j in c("te","tm")){
      for(i in c("02","04","08","16","32","64")){
         if(k=="bbp"){
            if(j=="tm"){
               fli<-paste("arc/",k,"/dat/",k,".vswf.mp.",i,".Rdata",sep="")
               flo<-paste("arc/",k,"/out/",k,".vswf.mp.",i,".dat"  ,sep="")
               cat(k,"mp",i,"\n")
               fli<-paste("arc/",k,"/dat/",k,".vswf.mp.",i,".Rdata",sep="")
               flo<-paste("arc/",k,"/out/",k,".vswf.mp.",i,".dat"  ,sep="")
            }else{
               fli<-paste("arc/",k,"/dat/",k,".vswf.pp.",i,".Rdata",sep="")
               flo<-paste("arc/",k,"/out/",k,".vswf.pp.",i,".dat"  ,sep="")
               cat(k,"pp",i,"\n")
               fli<-paste("arc/",k,"/dat/",k,".vswf.pp.",i,".Rdata",sep="")
               flo<-paste("arc/",k,"/out/",k,".vswf.pp.",i,".dat"  ,sep="")
            }
         }else{
            fli<-paste("arc/",k,"/dat/",k,".vswf.",j,".",i,".Rdata",sep="")
            flo<-paste("arc/",k,"/out/",k,".vswf.",j,".",i,".dat"  ,sep="")
            cat(k,j,i,"\n")
            fli<-paste("arc/",k,"/dat/",k,".vswf.",j,".",i,".Rdata",sep="")
            flo<-paste("arc/",k,"/out/",k,".vswf.",j,".",i,".dat"  ,sep="")
         }
         if(file.exists(fli)){
            load(fli)
            df<-data.frame(x=PWE$rx,
                           y=PWE$ry,
                           z=PWE$rz,
                           ReHm=Re(PWE$Hm),
                           ReHz=Re(PWE$Hz),
                           ReHp=Re(PWE$Hp),
                           ReEm=Re(PWE$Em),
                           ReEz=Re(PWE$Ez),
                           ReEp=Re(PWE$Ep),
                           ImHm=Im(PWE$Hm),
                           ImHz=Im(PWE$Hz),
                           ImHp=Im(PWE$Hp),
                           ImEm=Im(PWE$Em),
                           ImEz=Im(PWE$Ez),
                           ImEp=Im(PWE$Ep))
            write.fwf(file=flo,df,rownames=FALSE)
         }else{
            cat(fli,"nao existe!\n")
         }
      }
      if(k=="bbp"){
         if(j=="tm"){
            fli<-paste("arc/",k,"/dat/",k,".vwfd.mp.00.Rdata",sep="")
            flo<-paste("arc/",k,"/out/",k,".vwfd.mp.00.dat"  ,sep="")
         }else{
            fli<-paste("arc/",k,"/dat/",k,".vwfd.pp.00.Rdata",sep="")
            flo<-paste("arc/",k,"/out/",k,".vwfd.pp.00.dat"  ,sep="")
         }
      }else{
         fli<-paste("arc/",k,"/dat/",k,".vwfd.",j,".00.Rdata",sep="")
         flo<-paste("arc/",k,"/out/",k,".vwfd.",j,".00.dat"  ,sep="")
      }
      load(fli)
      if(k=="cwg"){DEF<-CWG}
      if(k=="rwg"){DEF<-RWG}
      if(k=="bbz"){DEF<-BBZ}
      if(k=="bbp"){DEF<-BBP}
      if(file.exists(fli)){
         load(fli)
         df<-data.frame(
                  x=DEF$rx,
                  y=DEF$ry,
                  z=DEF$rz,
                  ReHm=Re(DEF$Hm),
                  ReHz=Re(DEF$Hz),
                  ReHp=Re(DEF$Hp),
                  ReEm=Re(DEF$Em),
                  ReEz=Re(DEF$Ez),
                  ReEp=Re(DEF$Ep),
                  ImHm=Im(DEF$Hm),
                  ImHz=Im(DEF$Hz),
                  ImHp=Im(DEF$Hp),
                  ImEm=Im(DEF$Em),
                  ImEz=Im(DEF$Ez),
                  ImEp=Im(DEF$Ep))
         write.fwf(file=flo,df,rownames=FALSE)
      }else{
         cat(fli,"nao existe!\n")
      }
   }
}
