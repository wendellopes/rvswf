#-------------------------------------------------------------------------------
# Calculo do vetor de Poynting VWFD
#-------------------------------------------------------------------------------
Sm.BBP<-1i*(BBP$Ez*Conj(BBP$Hp)-BBP$Ep*Conj(BBP$Hz))/2
Sz.BBP<-1i*(BBP$Ep*Conj(BBP$Hp)-BBP$Em*Conj(BBP$Hm))/2
Sp.BBP<-1i*(BBP$Em*Conj(BBP$Hz)-BBP$Ez*Conj(BBP$Hm))/2
#-------------------------------------------------------------------------------
x11();image(BBP$k*BBP$y/(2*pi),BBP$k*BBP$x/(2*pi),abs(array(Sz.BBP,c(BBP$ny,BBP$nx))),main="abs Sz",col=cm.colors(1024))
#    contour(BBP$k*BBP$x/(2*pi),BBP$k*BBP$y/(2*pi),abs(array(Sz.BBP,c(BBP$nx,BBP$ny))),add=TRUE,drawlabels=FALSE)
x11();image(BBP$k*BBP$y/(2*pi),BBP$k*BBP$x/(2*pi),Re(array(Sz.BBP,c(BBP$ny,BBP$nx))),main=" Re Sz",col=cm.colors(1024))
x11();image(BBP$k*BBP$y/(2*pi),BBP$k*BBP$x/(2*pi),Im(array(Sz.BBP,c(BBP$ny,BBP$nx))),main=" Im Sz",col=cm.colors(1024))
#-------------------------------------------------------------------------------
# Calculo do vetor de Poynting VSWF
#-------------------------------------------------------------------------------
Sm.PWE<-1i*(PWE$Ez*Conj(PWE$Hp)-PWE$Ep*Conj(PWE$Hz))/2
Sz.PWE<-1i*(PWE$Ep*Conj(PWE$Hp)-PWE$Em*Conj(PWE$Hm))/2
Sp.PWE<-1i*(PWE$Em*Conj(PWE$Hz)-PWE$Ez*Conj(PWE$Hm))/2
#-------------------------------------------------------------------------------
x11();image(PWE$k*PWE$y/(2*pi),PWE$k*PWE$x/(2*pi),abs(array(Sz.PWE,c(PWE$ny,PWE$nx))),main="abs Sz",col=cm.colors(1024))
#    contour(PWE$k*PWE$x/(2*pi),PWE$k*PWE$y/(2*pi),abs(array(Sz.PWE,c(PWE$nx,PWE$ny))),add=TRUE,drawlabels=FALSE)
x11();image(PWE$k*PWE$y/(2*pi),PWE$k*PWE$x/(2*pi),Re(array(Sz.PWE,c(PWE$ny,PWE$nx))),main=" Re Sz",col=cm.colors(1024))
x11();image(PWE$k*PWE$y/(2*pi),PWE$k*PWE$x/(2*pi),Im(array(Sz.PWE,c(PWE$ny,PWE$nx))),main=" Im Sz",col=cm.colors(1024))

