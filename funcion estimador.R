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

rmax <- 0.2
P<-closepairs(X,0.2,twice = TRUE)
#Supuestos:
#Homogeneidad
#Estacioneriedad
#las funcion jac pide el proceso puntual (para calcular la sumatoria) y un conjunto de puntos (para calcular la integral por el metodo montecarlo)
#la funcion jac devuelve una matriz 
#la funciones pcf, y est devuelven un vector
#la funcion est pide la cordenada que se desea calcular
#Mejora entregar dos listas de parametros uno que sean lineales y otros no lineales

lavn<- function(X,pcf,est,jac,inipar,tol,kmax,mc=1000){
  x0<-runifpoint(mc,win=W)
  y0<-runifpoint(mc,win=W)
  dU<-crosspairs(x0,y0,0.2,what="ijd")
  P<-closepairs(X,0.2,twice = TRUE)
  par=inipar
  for(i in 1:kmax){
    r=c()
    for (i in 1:length(par)) {
      s=sum(dpcf(P$d,par,i))
      int=1/mc*1/mc*sum(pcf(dU$d,par)*est(dU$d,par,i))*intensity(X)^2
      e=s-int
      r=c(r,e)
    }
    gra=r
    jacob=jac(X,P$d,dU$d,par)
    jacob=solve(jacob)
    para=par
    par=par-(jacob%*%gra)[,1]
    e=sqrt(sum((para-par)*(para-par)))/sqrt(sum(para*para))
    if(e<tol){
      return(par)}
    print(par)
  }
  return(par)}

pcfe<- function(d,par){
  e=exp(par[1]^2*exp(-d/exp(par[2])))
  return(e)
}

jaco<-function(X,dp,du,par){
  cove<- function(d,par2){
    e=exp(-d/exp(par2))
    return(e)
  }
  rho=intensity(X)^2
  par1=par[1]
  par2=par[2]
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
j21=j11
j22=par1^2*s3-rho*par1^2*i22
jac=matrix(c(j11,j21,j12,j22),nrow=2,ncol=2)
return(jac)}

dpcf<- function(d,par,j){
    r=0
   if(j==1){
     r=2*par[1]*exp(-d/exp(par[2]))
   }  
   else if (j==2){
     r=2*par[1]*exp(-d/exp(par[2]))*exp(-par[2])*d
   }
  return(r)
}


mc=10
x0<-runifpoint(mc,win=W)
y0<-runifpoint(mc,win=W)
dU<-crosspairs(x0,y0,0.2,what="ijd")
dU$d

par=c(1,0.1)
pcfe(dU$d,par)
lavn(X,pcfe,dpcf,jaco,c(0.4,log(0.5)),tol=0.000000000000001,200,mc=1000)
log(0.15)
