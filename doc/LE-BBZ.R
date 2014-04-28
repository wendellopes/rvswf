#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
#source("BesselBeamsZ.r")
#source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
BBZ<-BesselBeamsZ(TE=TRUE,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBZ<-array(BBZ$Em,c(NP,NP,1))[,,1]
tez.BBZ<-array(BBZ$Ez,c(NP,NP,1))[,,1]
tep.BBZ<-array(BBZ$Ep,c(NP,NP,1))[,,1]
thm.BBZ<-array(BBZ$Hm,c(NP,NP,1))[,,1]
thz.BBZ<-array(BBZ$Hz,c(NP,NP,1))[,,1]
thp.BBZ<-array(BBZ$Hp,c(NP,NP,1))[,,1]
#
pdf("ReTemBbzWfd.pdf");image(x+xo,y+yo,Re(tem.BBZ),main="BBZ Em",col=cm.colors(1024));grid();dev.off()
pdf("ReTezBbzWfd.pdf");image(x+xo,y+yo,Re(tez.BBZ),main="BBZ Ez",col=cm.colors(1024));grid();dev.off()
pdf("ReTepBbzWfd.pdf");image(x+xo,y+yo,Re(tep.BBZ),main="BBZ Ep",col=cm.colors(1024));grid();dev.off()
pdf("ReThmBbzWfd.pdf");image(x+xo,y+yo,Re(thm.BBZ),main="BBZ Hm",col=cm.colors(1024));grid();dev.off()
pdf("ReThzBbzWfd.pdf");image(x+xo,y+yo,Re(thz.BBZ),main="BBZ Hz",col=cm.colors(1024));grid();dev.off()
pdf("ReThpBbzWfd.pdf");image(x+xo,y+yo,Re(thp.BBZ),main="BBZ Hp",col=cm.colors(1024));grid();dev.off()
#
pdf("ImTemBbzWfd.pdf");image(x+xo,y+yo,Im(tem.BBZ),main="BBZ Em",col=cm.colors(1024));grid();dev.off()
pdf("ImTezBbzWfd.pdf");image(x+xo,y+yo,Im(tez.BBZ),main="BBZ Ez",col=cm.colors(1024));grid();dev.off()
pdf("ImTepBbzWfd.pdf");image(x+xo,y+yo,Im(tep.BBZ),main="BBZ Ep",col=cm.colors(1024));grid();dev.off()
pdf("ImThmBbzWfd.pdf");image(x+xo,y+yo,Im(thm.BBZ),main="BBZ Hm",col=cm.colors(1024));grid();dev.off()
pdf("ImThzBbzWfd.pdf");image(x+xo,y+yo,Im(thz.BBZ),main="BBZ Hz",col=cm.colors(1024));grid();dev.off()
pdf("ImThpBbzWfd.pdf");image(x+xo,y+yo,Im(thp.BBZ),main="BBZ Hp",col=cm.colors(1024));grid();dev.off()

par(mfrow(6,2))
image(x+xo,y+yo,Re(tem.BBZ),main="BBZ Em",col=cm.colors(1024));grid()
image(x+xo,y+yo,Re(tez.BBZ),main="BBZ Ez",col=cm.colors(1024));grid()
image(x+xo,y+yo,Re(tep.BBZ),main="BBZ Ep",col=cm.colors(1024));grid()
image(x+xo,y+yo,Re(thm.BBZ),main="BBZ Hm",col=cm.colors(1024));grid()
image(x+xo,y+yo,Re(thz.BBZ),main="BBZ Hz",col=cm.colors(1024));grid()
image(x+xo,y+yo,Re(thp.BBZ),main="BBZ Hp",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(tem.BBZ),main="BBZ Em",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(tez.BBZ),main="BBZ Ez",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(tep.BBZ),main="BBZ Ep",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(thm.BBZ),main="BBZ Hm",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(thz.BBZ),main="BBZ Hz",col=cm.colors(1024));grid()
image(x+xo,y+yo,Im(thp.BBZ),main="BBZ Hp",col=cm.colors(1024));grid()

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
pdf("ReTemBbzPwe.pdf");image(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTezBbzPwe.pdf");image(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTepBbzPwe.pdf");image(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThmBbzPwe.pdf");image(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThzBbzPwe.pdf");image(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThpBbzPwe.pdf");image(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
# Parte Imaginaria
pdf("ImTemBbzPwe.pdf");image(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTezBbzPwe.pdf");image(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTepBbzPwe.pdf");image(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThmBbzPwe.pdf");image(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThzBbzPwe.pdf");image(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThpBbzPwe.pdf");image(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------



