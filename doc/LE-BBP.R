#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
#source("BesselBeamsP.r")
#source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
BBP<-BesselBeamsP(P,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBP<-array(BBP$Em,c(NP,NP,1))[,,1]
tez.BBP<-array(BBP$Ez,c(NP,NP,1))[,,1]
tep.BBP<-array(BBP$Ep,c(NP,NP,1))[,,1]
thm.BBP<-array(BBP$Hm,c(NP,NP,1))[,,1]
thz.BBP<-array(BBP$Hz,c(NP,NP,1))[,,1]
thp.BBP<-array(BBP$Hp,c(NP,NP,1))[,,1]
# Parte Real
pdf("ReTemBbpWfd.pdf");image(x+xo,y+yo,Re(tem.BBP),main="BBP Em",col=cm.colors(1024));grid();dev.off()
pdf("ReTezBbpWfd.pdf");image(x+xo,y+yo,Re(tez.BBP),main="BBP Ez",col=cm.colors(1024));grid();dev.off()
pdf("ReTepBbpWfd.pdf");image(x+xo,y+yo,Re(tep.BBP),main="BBP Ep",col=cm.colors(1024));grid();dev.off()
pdf("ReThmBbpWfd.pdf");image(x+xo,y+yo,Re(thm.BBP),main="BBP Hm",col=cm.colors(1024));grid();dev.off()
pdf("ReThzBbpWfd.pdf");image(x+xo,y+yo,Re(thz.BBP),main="BBP Hz",col=cm.colors(1024));grid();dev.off()
pdf("ReThpBbpWfd.pdf");image(x+xo,y+yo,Re(thp.BBP),main="BBP Hp",col=cm.colors(1024));grid();dev.off()
# ParteTmaBnaWa
pdf("ImTemBbpWfd.pdf");image(x+xo,y+yo,Im(tem.BBP),main="BBP Em",col=cm.colors(1024));grid();dev.off()
pdf("ImTezBbpWfd.pdf");image(x+xo,y+yo,Im(tez.BBP),main="BBP Ez",col=cm.colors(1024));grid();dev.off()
pdf("ImTepBbpWfd.pdf");image(x+xo,y+yo,Im(tep.BBP),main="BBP Ep",col=cm.colors(1024));grid();dev.off()
pdf("ImThmBbpWfd.pdf");image(x+xo,y+yo,Im(thm.BBP),main="BBP Hm",col=cm.colors(1024));grid();dev.off()
pdf("ImThzBbpWfd.pdf");image(x+xo,y+yo,Im(thz.BBP),main="BBP Hz",col=cm.colors(1024));grid();dev.off()
pdf("ImThpBbpWfd.pdf");image(x+xo,y+yo,Im(thp.BBP),main="BBP Hp",col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
UVS<-PVectorSphericalWaveFunctions(k,x,y,ze,lmax,u$GTE,u$GTM)
#
tem.uvs<-array(UVS$Em,c(NP,NP,1))[,,1]
tez.uvs<-array(UVS$Ez,c(NP,NP,1))[,,1]
tep.uvs<-array(UVS$Ep,c(NP,NP,1))[,,1]
thm.uvs<-array(UVS$Hm,c(NP,NP,1))[,,1]
thz.uvs<-array(UVS$Hz,c(NP,NP,1))[,,1]
thp.uvs<-array(UVS$Hp,c(NP,NP,1))[,,1]
# Parte Real
pdf("ReTemBbpPwe.pdf");image(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTezBbpPwe.pdf");image(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTepBbpPwe.pdf");image(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThmBbpPwe.pdf");image(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThzBbpPwe.pdf");image(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThpBbpPwe.pdf");image(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
# ParteTmaBnaWa
pdf("ImTemBbpPwe.pdf");image(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTezBbpPwe.pdf");image(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTepBbpPwe.pdf");image(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThmBbpPwe.pdf");image(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThzBbpPwe.pdf");image(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThpBbpPwe.pdf");image(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------
