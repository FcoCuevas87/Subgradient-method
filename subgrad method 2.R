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
#se calcula el close pairs antes
P<-closepairs(X,rmax=1,twice = FALSE)
#integral en el estimador
inte<-function(arg,par,par2,cova,fgrad,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cova(d,par,par2))*fgrad(d,par,par2,j)
  return(w)
}
#gradiente
est<-function(par,par2,f,cova,j){
  int=vegas(inte, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,par2=par2,cova=cova,fgrad=f,j=j)
  a=f(P$d,par,par2,j)
  a=sum(a)-int$int
  return(-a)
}
#Waageralgo
fwag<-function(d,par,par2,j){
  return(cos(par2[j]*d))*(1*(d<=0.2) + 0* (d>=0.2))
}
cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
    a=a+par[i]*cos(par2[i]*d)
    # falta la funcion k
  }
  return(a)
}
par=c(0.0005,0.0005,0.0005,0.0005,0.0005)
par2=c(0.05,0.1,0.15,0.2,0.25)
alp=10
#k2=0.5
lam=0.25
est(para,par2,fwag,cove,1)
for (k in 1:100){
  para=par
  print(par)
  for (i in 1:(length(par))){
    gr=-est(para,par2,fwag,cove,i)
    d=para[i]
    if (abs(para[i]-gr)>lam){
      d=-gr+sign(para[i]-gr)*lam
    }
    par[i]=para[i]+alp*d
  }
}

