rm(list = ls())
library(spatstat)
library(cubature)
#simulacion LGCP
N=500
sig=1
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)

Z<-attr(X,"Lambda")
plot(Z)
points(X,pch=16)
X
rmax <- 0.2
P<-closepairs(X,0.2,twice = TRUE)

#funcion de covarianza util para simpolificar el calculo
cove<- function(d,par2){
  e=exp(-d/exp(par2))
  return(e)
}
#integrales que aparecen en sus respectivos terminos
inte1<- function(d,par1,par2){
  re=cove(d,par2)*exp(par1^2*cove(d,par2))
  re=sum(re)
  return(re)
}
inte2<-function(d,par1,par2){
  re=exp(par1^2*cove(d,par2))*cove(d,par2)*exp(-par2)*d
  re=sum(re)
  return(re)
}

inte11<-function(d,par1,par2){
  re=exp(par1^2*cove(d,par2))*cove(d,par2)*(1+2*par1^2*cove(d,par2))
  re=sum(re)
  return(re)
}

inte12<-function (d,par1,par2){
  re=exp(par1^2*cove(d,par2))*cove(d,par2)*exp(-par2)*d*(1+par1^2*cove(d,par2))
  re=sum(re)
  return(re)
}
inte21<-function(d,par1,par2){
  re=exp(par1^2*cove(d,par2))*cove(d,par2)*exp(-par2)*d*(1+par1^2*cove(d,par2))
  re=sum(re)
  return(re)
}
inte22<-function(d,par1,par2){
  re=exp(par1^2*cove(d,par2))*cove(d,par2)*exp(-par2)*d*(cove(d,par2)*exp(-par2)*d*par1^2+exp(-par2)*d-1)
  re=sum(re)
  return(re)
}
#Sumas que se repiten en varios terminos

sum1<-function(d,par2){
  re=cove(d,par2)
  re=sum(re)
  return(re)
}

sum2<-function(d,par2){
  re=cove(d,par2)*exp(-par2)*d
  re=sum(re)
  return(re)
}

sum3<-function(d,par2){
  re=cove(d,par2)*exp(-par2)*d*(exp(-par2)*d-1)
  re=sum(re)
  return(re)
}


mc=1000
rho=intensity(X)^2
x0<-runifpoint(mc,win=W)
y0<-runifpoint(mc,win=W)
dU<-crosspairs(x0,y0,rmax,what="ijd")
P<-closepairs(X,0.2,twice = TRUE)
par1=1
par2=log(0.15)

#Newton Rapson varianza
for (i in 1:200){
  #obtenemos los valores de las integrales
  i1=(1/mc)^2*inte1(dU$d,par1,par2)
  i2=(1/mc)^2*inte2(dU$d,par1,par2)
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i21=(1/mc)^2*inte21(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos los valores de los gradientes
  e1=2*par1*s1-2*rho*par1*i1
  e2=par1^2*s2-rho*par1^2*i2
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho*i11
  par1=par1-e1/j11
  print(par1^2)
}

#Newton Rapson parametro2
for (i in 1:200){
  #obtenemos los valores de las integrales
  i1=(1/mc)^2*inte1(dU$d,par1,par2)
  i2=(1/mc)^2*inte2(dU$d,par1,par2)
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i21=(1/mc)^2*inte21(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos los valores de los gradiente
  e2=par1^2*s2-rho*par1^2*i2
  #remplazamos en los valores de la jacobiana
  j22=par1^2*s3-rho*par1^2*i22
  par2=par2-e2/j22
  print(c(exp(par2)))
}


mc=1000
rho=intensity(X)^2
x0<-runifpoint(mc,win=W)
y0<-runifpoint(mc,win=W)
dU<-crosspairs(x0,y0,rmax,what="ijd")
P<-closepairs(X,0.2,twice = TRUE)
par1=1
par2=log(0.15)

#Newton Rapson Jacobiana

for (i in 1:200){
  #obtenemos los valores de las integrales
  i1=(1/mc)^2*inte1(dU$d,par1,par2)
  i2=(1/mc)^2*inte2(dU$d,par1,par2)
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i21=(1/mc)^2*inte21(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos los valores de los gradientes
  e1=2*par1*s1-2*rho*par1*i1
  e2=par1^2*s2-rho*par1^2*i2
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho*i11
  j12=2*par1*s2-2*rho*par1*i12
  j21=2*par1*s2-2*rho*par1*i21
  j22=par1^2*s3-rho*par1^2*i22
  c(j11,j21,j12,j22)
  jac=matrix(c(j11,j21,j12,j22),nrow=2,ncol=2)
  jac=solve(jac)
  x=c(par1,par2)
  x=x-(jac%*%c(e1,e2))[,1]
  par1=x[1]
  par2=x[2]
  print(c(par1^2,exp(par2)))
}

#Newton Rapson por cordenadas (o usando la diagonal)
for (i in 1:1000){
  #obtenemos los valores de las integrales
  i1=(1/mc)^2*inte1(dU$d,par1,par2)
  i2=(1/mc)^2*inte2(dU$d,par1,par2)
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i21=(1/mc)^2*inte21(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos los valores de los gradientes
  e1=2*par1*s1-rho*par1*i1
  e2=par1^2*s2-rho*par1^2*i2
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho^2*i11
  j22=par1^2*s3-rho*par1^2*i22
  c(j11,j21,j12,j22)
  jac=matrix(c(j11,0,0,j22),nrow=2,ncol=2)
  jac
  jac=solve(jac)
  x=c(par1,par2)
  x=x-(jac%*%c(e1,e2))[,1]
  par1=x[1]
  par2=x[2]
  print(c(par1^2,exp(par2)))
}


c(par1^2,exp(par2))
jac=matrix(c(11,21,12,22),nrow=2,ncol=2)
jac