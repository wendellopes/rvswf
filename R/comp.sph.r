x<-5
nmax<-10
n<-0:10
v.cal<-bess.sph(x,nmax)$jn
v.ref<-reff.sjn(x,0:nmax)
print(cbind(n,v.ref,v.cal))
