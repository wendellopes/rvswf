#-------------------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
# Imagem com legenda
#-------------------------------------------------------------------------------
source("~/R-VSWF/opex/Image.R")
#-------------------------------------------------------------------------------
# Carregando os valores
#-------------------------------------------------------------------------------
#load("rwg.vwfd.tm.00.Rdata")
load("rwg.vswf.tm.16.Rdata")
print(ls())
RWG<-PWE
#-------------------------------------------------------------------------------
# Primeira olhada
#-------------------------------------------------------------------------------
#str(RWG)
#-------------------------------------------------------------------------------
# Calculo do vetor de Poynting
#-------------------------------------------------------------------------------
Sm<-1i*(RWG$Ez*Conj(RWG$Hp)-RWG$Ep*Conj(RWG$Hz))/2
Sz<-1i*(RWG$Ep*Conj(RWG$Hp)-RWG$Em*Conj(RWG$Hm))/2
Sp<-1i*(RWG$Em*Conj(RWG$Hz)-RWG$Ez*Conj(RWG$Hm))/2
#-------------------------------------------------------------------------------
# USANDO X e Y 
#-------------------------------------------------------------------------------
# Sx<-(Sp+Sm)/(   sqrt(2))
# Sy<-(Sp-Sm)/(1i*sqrt(2))
# Sz<-Sz
# # Campos em xy
# Ex<-(RWG$Ep+RWG$Em)/(   sqrt(2))
# Ey<-(RWG$Ep-RWG$Em)/(1i*sqrt(2))
# Ez<-RWG$Ez
# 
# Hx<-(RWG$Hp+RWG$Hm)/(   sqrt(2))
# Hy<-(RWG$Hp-RWG$Hm)/(1i*sqrt(2))
# Hz<-RWG$Hz
# # Vetor de Poynting
# SSx<-(Ey*Conj(Hz)-Ez*Conj(Hy))/2
# SSy<-(Ez*Conj(Hx)-Ex*Conj(Hz))/2
# SSz<-(Ex*Conj(Hy)-Ey*Conj(Hx))/2
#-------------------------------------------------------------------------------
# Comparando
#-------------------------------------------------------------------------------
## Partes reais
#range(Re( Sx))
#range(Re(SSx))
#range(Re( Sy))
#range(Re(SSy))
#range(Re( Sz))
#range(Re(SSz))
## Partes Imaginarias
#range(Im( Sx))
#range(Im(SSx))
#range(Im( Sy))
#range(Im(SSy))
#range(Im( Sz))
#range(Im(SSz))
#-------------------------------------------------------------------------------
# GRAFICOS
#-------------------------------------------------------------------------------
# # partes reais e imaginarias
# x11();Image(Re(array(Sz,c(200,200))),main="Re Sz",,col=cm.colors(1024))
# x11();Image(Im(array(Sz,c(200,200))),main="Im Sz",,col=cm.colors(1024))
# x11();Image(Re(array(Sp,c(200,200))),main="Re Sp",,col=cm.colors(1024))
# x11();Image(Im(array(Sp,c(200,200))),main="Im Sp",,col=cm.colors(1024))
# x11();Image(Re(array(Sm,c(200,200))),main="Re Sm",,col=cm.colors(1024))
# x11();Image(Im(array(Sm,c(200,200))),main="Im Sm",,col=cm.colors(1024))
# # Amplitudes
x11();image(RWG$k*RWG$x/(2*pi),RWG$k*RWG$y/(2*pi),abs(array(Sz,c(200,200))),main="abs Sz",col=heat.colors(1024))
contour(RWG$k*RWG$x/(2*pi),RWG$k*RWG$y/(2*pi),abs(array(Sz,c(200,200))),add=TRUE,drawlabels=FALSE)
# x11();image(abs(array(Sx,c(200,200))),main="abs Sx",col=cm.colors(1024))
# x11();image(abs(array(Sy,c(200,200))),main="abs Sy",col=cm.colors(1024))
# #
# x11();image(abs(array(SSz,c(200,200))),main="abs SSz",col=cm.colors(1024))
# x11();image(abs(array(SSx,c(200,200))),main="abs SSx",col=cm.colors(1024))
# x11();image(abs(array(SSy,c(200,200))),main="abs SSy",col=cm.colors(1024))


