library(spatstat)
library(cubature)
#simulacion LGCP
N=500
sig=0.1
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
rmax=0.2
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)

Z<-attr(X,"Lambda")
plot(Z)
points(X,pch=16)
X

P<-closepairs(X,0.2,twice = FALSE)
cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
    a=a+par[i]*exp(-d/par2[i]) 
    # falta la funcion k
  }
  return(a)
}
fwag<- function(d,par,par2,j){
  return(par[j]*exp(-d/par2[j])*(1*(d<=0.2) + 0* (d>=0.2)))
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
cont=1
par=runif(5, 0, 1)
par=c(0.5,0.5,0.5,0.5,0.5)
par2=c(0.05,0.1,0.15,0.2,0.25)
mc=1000
#no hago un for pq sino mi ram explota
sum=0
result=vegas(integrand1, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
                   relTol=1e-3, absTol=1e-6,
                   flags=list(verbose=0, final=1),par=par,par2=par2, j=3)
for (i in 1:20){
x0<-runifpoint(mc,win=square(1))
y0<-runifpoint(mc,win=square(1))
dU<-crosspairs(x0,y0,rmax,what="ijd")
resmc=(1/mc)^2*sum(integrandM(dU$d,par,par2,3))
dif=abs(result$int-resmc)/abs(result$int)*100
sum=sum+dif
print(dif)}
print(sum/40)

