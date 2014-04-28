#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
#source("CylindricalWaveGuide.r")
#source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
CWG<-CylindricalWaveGuide(TE=FALSE,M,S,g,kz,x+xo,y+yo,ze+zo)
#
tem.cwg<-array(CWG$Em,c(NP,NP,1))[,,1]
tez.cwg<-array(CWG$Ez,c(NP,NP,1))[,,1]
tep.cwg<-array(CWG$Ep,c(NP,NP,1))[,,1]
thm.cwg<-array(CWG$Hm,c(NP,NP,1))[,,1]
thz.cwg<-array(CWG$Hz,c(NP,NP,1))[,,1]
thp.cwg<-array(CWG$Hp,c(NP,NP,1))[,,1]
# Parte Real
pdf("ReTemCwgWfd.pdf");image(x+xo,y+yo,Re(tem.cwg),main="CWG Em",col=cm.colors(1024));grid();dev.off()
pdf("ReTezCwgWfd.pdf");image(x+xo,y+yo,Re(tez.cwg),main="CWG Ez",col=cm.colors(1024));grid();dev.off()
pdf("ReTepCwgWfd.pdf");image(x+xo,y+yo,Re(tep.cwg),main="CWG Ep",col=cm.colors(1024));grid();dev.off()
pdf("ReThmCwgWfd.pdf");image(x+xo,y+yo,Re(thm.cwg),main="CWG Hm",col=cm.colors(1024));grid();dev.off()
pdf("ReThzCwgWfd.pdf");image(x+xo,y+yo,Re(thz.cwg),main="CWG Hz",col=cm.colors(1024));grid();dev.off()
pdf("ReThpCwgWfd.pdf");image(x+xo,y+yo,Re(thp.cwg),main="CWG Hp",col=cm.colors(1024));grid();dev.off()
# Parte Imaginaria
pdf("ImTemCwgWfd.pdf");image(x+xo,y+yo,Im(tem.cwg),main="CWG Em",col=cm.colors(1024));grid();dev.off()
pdf("ImTezCwgWfd.pdf");image(x+xo,y+yo,Im(tez.cwg),main="CWG Ez",col=cm.colors(1024));grid();dev.off()
pdf("ImTepCwgWfd.pdf");image(x+xo,y+yo,Im(tep.cwg),main="CWG Ep",col=cm.colors(1024));grid();dev.off()
pdf("ImThmCwgWfd.pdf");image(x+xo,y+yo,Im(thm.cwg),main="CWG Hm",col=cm.colors(1024));grid();dev.off()
pdf("ImThzCwgWfd.pdf");image(x+xo,y+yo,Im(thz.cwg),main="CWG Hz",col=cm.colors(1024));grid();dev.off()
pdf("ImThpCwgWfd.pdf");image(x+xo,y+yo,Im(thp.cwg),main="CWG Hp",col=cm.colors(1024));grid();dev.off()
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
pdf("ReTemCwgPwe.pdf");image(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTezCwgPwe.pdf");image(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReTepCwgPwe.pdf");image(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThmCwgPwe.pdf");image(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThzCwgPwe.pdf");image(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ReThpCwgPwe.pdf");image(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
# Parte Imaginaria
pdf("ImTemCwgPwe.pdf");image(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTezCwgPwe.pdf");image(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImTepCwgPwe.pdf");image(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThmCwgPwe.pdf");image(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThzCwgPwe.pdf");image(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf("ImThpCwgPwe.pdf");image(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------
