rm(list = ls())

x=matrix(rnorm(1000*5,0,1),ncol=5,nrow=1000)
x=cbind(1,x)
b=c(2,7,8,12,1,-2)
y=x%*%b+rnorm(1000,0,1)
lm.fit<-lm(y~x+0)
lm.fit
#


grad<-function(x,y,b,i){
  return(t(x[,i])%*%(x%*%b-y))
}

grad(x,y,b,4)
t(x[,i])%*%(x%*%b-y)
par=c(0.1,0.2,0.3,0.4,0.5,0.6)
lam=10000
alp=0.001
for (k in 1:100){
  para=par
  print(par)
  for (i in 1:(length(par))){
    x0=para[i]-alp*grad(x,y,par,i)
    par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)
  }
}


