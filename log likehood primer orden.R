library(spatstat)
library(cubature)
inten<- function(x,y) (25*x^2+16*x+4*y+9*y^2)
X <-rpoispp(100,lambda=inten)
plot(X)
integrand<-function(arg,j){
  x<-arg[1]
  y<-arg[2]
  d=c(0,0,0,0,0)
  d[j]=1
  w=d[1]*x^2+d[2]*x^1+d[3]+d[4]*y+d[5]*y^2
  return(w)
}
grad<-function(X,b,j){
  int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
            relTol=1e-3, absTol=1e-12,
            flags=list(verbose=0, final=1) ,j=j)
return(X$n/b[j]-int$int)
}
int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
          relTol=1e-3, absTol=1e-12,
          flags=list(verbose=0, final=1) ,j=1)
int$int
par=c(sqrt(10),sqrt(10),sqrt(10),sqrt(10),sqrt(10))
lam=0
alp=0.1
grad(X,par,5)
for (k in 1:200){
  para=par
  print(par^2)
  for (i in 1:(length(par))){
    x0=para[i]-alp*grad(X,par,i)
    par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)
  }
}

#resultado esperado
result=c(0,0,0,0,0)
for (i in 1:(length(result))){
  int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
            relTol=1e-3, absTol=1e-12,
            flags=list(verbose=0, final=1) ,j=i)
  result[i]=X$n/int$int
}
result
gplt<-function(b){
  int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
            relTol=1e-3, absTol=1e-12,
            flags=list(verbose=0, final=1) ,j=3)
  print(int$int)
  return(X$n/b-int$int)
}
