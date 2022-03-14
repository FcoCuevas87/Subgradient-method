library(spatstat)
library(cubature)
N=200
sig=0.5
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)
Z<-attr(X,"Lambda")
#covarianza=varianza^2*exp(-d*escala)
plot(Z)
points(X,pch=16)
P<-closepairs(X,rmax=1,twice = FALSE)
cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
  a=a+par[i]*cos(par2[i]*d) 
  # falta la funcion k
  }
  return(a)
}

integrand1<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

integrand2<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*cos(par2[j]*d)*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

grad<-function(par,par2,j,rmax){
  int=vegas(integrand1, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,par2=par2,j=j)
  int2=vegas(integrand2, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
             relTol=1e-3, absTol=1e-6,
             flags=list(verbose=0, final=1),par=par,par2=par2,j=j)
  P<-closepairs(X,rmax,twice = FALSE)
  a=cos(par2[j]*P$d)
  ver=a-int$int/int2$int
  ver=-sum(ver)
  return(ver)
}
par=c(1,2,3,4,5,6,7)
par2=c(2,3,4,5,6,7,8)
grad(par,par2,2,1)

#para cualquier f

inte<-function(arg,par,par2,cova,fgrad,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cova(d,par,par2))*fgrad(d,par,par2,j)
  return(w)
}
#Waageralgo
fwag<-function(d,par,par2,j){
  return(cos(par2[j]*d)*(1*(d<=0.2) + 0* (d>=0.2)))
    
}

#no se agrega el patron puntual para ahorrar tiempo
est<-function(par,par2,f,cova,j){
  int=vegas(inte, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,par2=par2,cova=cova,fgrad=f,j=j)
  a=f(P$d,par,par2,j)
  a=sum(a)-int$int
  return(-a)        
}



int=vegas(inte, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
          relTol=1e-3, absTol=1e-6,
          flags=list(verbose=0, final=1),par=par,par2=par2,cova=cove,fgrad=fwag,j=2)
a=fwag(P$d,par,par2,2)
sum(a)
int$int
est(par,par2,fwag,cove,2)
a=cos2(P$d,par,par2)
a=sum(a)-int$int
sum(a)

par=c(0.0005,0.0005,0.0005,0.0005,0.0005)
par2=c(0.05,0.1,0.15,0.2,0.25)
alp=0.5
k2=0.5
for (k in 1:20){
  para=par
  t=2
  print(par)
  for (i in 1:(length(par))){
     x=-est(para,par2,fwag,cove,i)+1/t*para[i]
     par[i]=t*sign(x)*max(c(abs(x)-alp*k2,0))/(t*(1+alp*(1-k2))+1)
  }
}

