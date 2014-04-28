\name{vswf.jlm}
\alias{vswf.jlm}
\title{Position j(l,m)}
\description{Gives the list position j(l,m)=l(l+1)+m+1. As we know, the
   first element in a vector in R is 1, so j(0,0)=1.}
\usage{
vswf.jlm(l,m)
}
\arguments{
  \item{l}{l>=0}
  \item{m}{-l<=m<=l}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
lmax=5
LMAX=vswf.jlm(lmax,lmax)
v<-rep(0,LMAX)
u<-data.frame(l=v,m=v,j=v)
for(l in 0:lmax){
   for(m in -l:l){
      j<-vswf.jlm(l,m)
      u$l[j]<-l
      u$m[j]<-m
      u$j[j]<-j
   }
}
print(u)
}
