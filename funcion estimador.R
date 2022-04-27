rm(list = ls())
library(spatstat)
library(cubature)
#simulacion LGCP
N=500
sig=3
W=square(1) ### 1.- Considerar el tamaño de la ventana en las integrales
Warea=spatstat.geom::area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W) ### (*).- Crear mi propio simulador de LGCP

Z<-attr(X,"Lambda")
plot(Z)
points(X,pch=16)

rmax <- 0.2 ### 2.- Añadirlo como parte de las funciones que necesiten closepairs y crosspairs
P<-closepairs(X,rmax,twice = TRUE) 
#Supuestos:
#Homogeneidad
#Estacioneriedad
#las funcion jac pide el proceso puntual (para calcular la sumatoria) y un conjunto de puntos (para calcular la integral por el metodo montecarlo)
#la funcion jac devuelve una matriz 
#la funciones pcf, y est devuelven un vector
#la funcion est pide la cordenada que se desea calcular
#Mejora entregar dos listas de parametros uno que sean lineales y otros no lineales

#Considerar diferentes w_{ij}
lavn<- function(X,g,est,jac,inipar,tol,kmax,mc=1000, trace = TRUE){
  x0<-runifpoint(mc,win=W)
  y0<-runifpoint(mc,win=W)
  dU<-crosspairs(x0,y0,rmax,what="ijd")
  P<-closepairs(X,rmax,twice = TRUE)
  par=inipar
  for(i in 1:kmax){
    r=c()
    for (j in 1:length(par)) {
      s=sum(est(P$d,par,j))
      int=1/mc*1/mc*sum(g(dU$d,par)*est(dU$d,par,j))*(intensity(X))^2
      e=s-int
      r=c(r,e)
    }
    gra=r
    jacob=jac(X,P$d,dU$d,par,mc)
    jacob=solve(jacob)
    para=par
    par=par-(jacob%*%gra)[,1]*0.1
    e=sqrt(sum((para-par)*(para-par)))/sqrt(sum(para*para))
    if(e<tol){
      return(par)}
    if(trace) cat(c(par[1]^2,exp(par[2])), i, "\n")
  }
  return(par)
}

pcfe<- function(d,par){
  e=exp(par[1]^2*exp(-d/exp(par[2])))
  return(e)
}

dpcf<- function(d,par,j){
  r=0
  if(j==1){
    r=2*par[1]*exp(-d/exp(par[2]))
  }  
  else if (j==2){
    r=par[1]^2*exp(-d/exp(par[2]))*exp(-par[2])*d
  }
  return(r)
}


jaco<-function(X,dp,du,par,mc){
  par1=par[1]
  par2=par[2]
  rho=(intensity(X))^2
  cove<- function(d,par2){
    e=exp(-d/exp(par2))
    return(e)
  }
  #integrales que aparecen en sus respectivos terminos
  
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
  i11=(1/mc)^2*inte11(du,par1,par2)
  i12=(1/mc)^2*inte12(du,par1,par2)
  i22=(1/mc)^2*inte22(du,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(dp,par2)
  s2=sum2(dp,par2)
  s3=sum3(dp,par2)
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho*i11
  j12=2*par1*s2-2*rho*par1*i12
  j21=j12
  j22=par1^2*s3-rho*par1^2*i22
  jac=matrix(c(j11,j21,j12,j22),nrow=2,ncol=2)
  return(jac)}


min.cv <- mincontrast(observed = pcf(X), theoretical = function(x,par) exp( (par[1]^2)*exp(-x/exp(par[2]))), startpar = c(sqrt(3),log(0.15)))
pcf(X)


help(mincontrast)
min.cv
s2 <- min.cv$par[1]
l.phi <- min.cv$par[2]
par=c(s2,l.phi)
re<-lavn(X,pcfe,dpcf,jaco,inipar=par,tol= 1e-15,kmax=200,mc=1000, trace = TRUE)




##### Considerar log(g) = \sum_{i = 1}^{n} sigma^2_{i}*exp(-h/phi_{i})

########### 2 testing ###########
#codigo para probar la función
cove<- function(d,par2){
  e=exp(-d/exp(par2))
  return(e)
}
#integrales que aparecen en sus respectivos terminos

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

kmax=200
mc=1000
x0<-runifpoint(mc,win=W)
y0<-runifpoint(mc,win=W)
dU<-crosspairs(x0,y0,rmax,what="ijd")
P<-closepairs(X,rmax,twice = TRUE)
par=c(0,log(0.15))
for(i in 1:kmax){
  r=c()
  for (i in 1:length(par)) {
    s=sum(dpcf(P$d,par,i))
    int=1/mc*1/mc*sum(pcfe(dU$d,par)*dpcf(dU$d,par,i))*(intensity(X))^2
    e=s-int
    r=c(r,e)
  }
  gra=r
  par1=par[1]
  par2=par[2]
  rho=(intensity(X))^2
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho*i11
  j12=2*par1*s2-2*rho*par1*i12
  j21=j12
  j22=par1^2*s3-rho*par1^2*i22
  jac2=jaco(X,P$d,dU$d,par,mc)
  jacob=matrix(c(j11,j21,j12,j22),nrow=2,ncol=2)
  jacob=solve(jacob)
  
  para=par
  par=par-(jacob%*%gra)[,1]
  e=sqrt(sum((para-par)*(para-par)))/sqrt(sum(para*para))
  print(c(par[1]^2,exp(par[2])))
  
}
print(gra)

par1=0.8
par2=log(0.15)

for (i in 1:200){
  #obtenemos los valores de las integrales
  
  i11=(1/mc)^2*inte11(dU$d,par1,par2)
  i12=(1/mc)^2*inte12(dU$d,par1,par2)
  i22=(1/mc)^2*inte22(dU$d,par1,par2)
  #obtenemos los valores de la sumatoria
  s1=sum1(P$d,par2)
  s2=sum2(P$d,par2)
  s3=sum3(P$d,par2)
  #remplazamos en los valores de la jacobiana
  j11=2*s1-2*rho*i11
  j12=2*par1*s2-2*rho*par1*i12
  j21=j12
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
