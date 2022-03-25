library(spatstat)
library(cubature)
#simulacion LGCP
N=500
sig=0.1
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)

Z<-attr(X,"Lambda")
plot(Z)
points(X,pch=16)
X

P<-closepairs(X,0.2,twice = FALSE)

#funciones necesarias

#Funcion de covarianza
cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
    a=a+exp(par[i])*exp(-d/par2[i]) 
    # falta la funcion k
  }
  return(a)
}
fwag<- function(d,par,par2,j){
  return(exp(par[i]*exp(-d/par2[i])*(1*(d<=0.2) + 0* (d>=0.2))))
}
#integral1 
integrand1<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*fwag(d,par,par2,j)
  return(w)
}
integrandM<-function(d,par,par2,j){
  w=exp(cove(d,par,par2))*fwag(d,par,par2,j)
  return(w)
}

#gradiente
grad<-function(par,par2,j){
  x0<-runifpoint(mc,win=square(1))
  y0<-runifpoint(mc,win=square(1))
  dU<-crosspairs(x0,y0,rmax,what="ijd")
  inte=(1/mc)^2*sum(integrandM(dU$d,par,par2,3))
  a=fwag(P$d,par,par2,j)
  ver=sum(a)-inte
  return(ver)
}
mc=1000
#varianzas (parametros a estimar)
par=c(0.5,0.5,0.5,0.5,0.5)
#phi
par2=c(0.05,0.1,0.15,0.2,0.25)
alp=0.000001
lam=0
for (k in 1:200){
  para=par
  print(par)
  for (i in 1:(length(par))){
    x0=para[i]-alp*grad(par,par2,i)
    par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)}
}

