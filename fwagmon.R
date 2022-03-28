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
  return(exp(par)*exp(-d/par2[j])*(1*(d<=0.2) + 0* (d>=0.2)))
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

integrandj<-function(d,par,par2,i,j){
  w=fwag(d,par,par2,i)*fwag(d,par,par2,j)*cove(d,par,par2)
  return(w)
}

#gradiente
gradien<-function(par,par2){
  g=vector()
  for (i in 1:length(par)){
    g[i]=grad(par,par2,i)
  }
  return(g)
}
grad<-function(par,par2,j){
  x0<-runifpoint(mc,win=square(1))
  y0<-runifpoint(mc,win=square(1))
  dU<-crosspairs(x0,y0,rmax,what="ijd")
  inte=(1/mc)^2*sum(integrandM(dU$d,par,par2,j))
  a=fwag(P$d,par,par2,j)
  ver=sum(a)-(intensity(X))^2*inte
  return(ver)
}
a=cbind(c(1,2),c(3,4))
b[1,]
jac=matrix(nrow=2,ncol=2)
jacob<-function(par,par2){
  x0<-runifpoint(mc,win=square(1))
  y0<-runifpoint(mc,win=square(1))
  dU<-crosspairs(x0,y0,rmax,what="ijd")
  jac=matrix(nrow=length(par),ncol=length(par))
  for (i in 1:length(par)){
    for (j in 1:length(par) ){
      temp=0
     temp=-(intensity(X))^2*(1/mc)^2*sum(integrandj(dU$d,par,par2,i,j))
     if(i==j){
       temp=temp+sum(fwag(P$d,par,par2,i))-(intensity(X))^2*(1/mc)^2*sum(integrandM(dU$d,par,par2,i))
     }
     jac[i,j]=temp
    }
  
  }
djacob<-function(par,par2){
  x0<-runifpoint(mc,win=square(1))
  y0<-runifpoint(mc,win=square(1))
  dU<-crosspairs(x0,y0,rmax,what="ijd")
  jac=c()
  for (i in 1:length(par)){
    temp=0
    temp=-(intensity(X))^2*(1/mc)^2*sum(integrandj(dU$d,par,par2,i,i))
    temp=temp+sum(fwag(P$d,par,par2,i))-(intensity(X))^2*(1/mc)^2*sum(integrandM(dU$d,par,par2,i))
    jac[i]=temp
    }
  return(jac)
}
c(1,2)*c(1,2)
1==2
mc=1000
#varianzas (parametros a estimar)
par=c(log(0.5),log(0.25))
#phi
par2=c(0.05,0.15)
alp=0.000001
lam=0
para=par
jac=jacob(par,par2)
jac
jac=solve(jac)
gr=gradien(par,par2)
gr
t(gr)
print(jac%*%gradien(par,par2))
for (k in 1:200){
  para=par
  jac=djacob(par,par2)
  jac=jac^{-1}
  x0=para-(jac*gradien(par,par2))
  par=x0
  print(exp(par))
}
a=c(2,4)
print(a^{-1})
for (k in 1:200){
  para=par
  jac=jacob(par,par2)
  jac=solve(jac)
  x0=para-(jac%*%gradien(par,par2))[,1]
  par=x0
  print(exp(par))
}
a=c(2,2)
a
ma=matrix(data=c(1,2,3,4),nrow=2,ncol=2)
ma
print(ma%*%a)
a=(ma%*%a)[,1]
a
for (k in 1:200){
  para=par
  print(par)
  for (i in 1:(length(par))){
    x0=para[i]-alp*grad(par,par2,i)
    par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)}
}

