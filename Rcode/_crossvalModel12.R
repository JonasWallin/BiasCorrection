##
# internal function for crossval.R
##

llike.mat2 <- function(p,obj)
{
  tau = exp(p[1])
  kappa = exp(p[2])
  mu = p[3]
  sigma2 = exp(p[4])

  if(tau>1e-16 && kappa>1e-16 && sigma2>1e-16)
  {
    Y_mu = obj$Y - mu
    Q0 = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q.p = kronecker(Diagonal(obj$n.cv-1,1),Q0)

    Q.p = obj$A%*%Q.p%*%obj$A
    AtQ <- t(obj$X)%*%Q.p
    Qhat <-  obj$I/sigma2  + AtQ%*%obj$X
    Rp <- chol(Qhat)
    v = solve(t(Rp),AtQ%*%Y_mu)
    l = (obj$n.cv-1)*sum(log(diag(chol(Q0))))-sum(log(diag(Rp))) - obj$N*log(sigma2)/2  + t(v)%*%v/2 - Y_mu%*%Q.p%*%Y_mu/2
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

llike.beta <- function(p,obj)
{
  tau = exp(p[1])
  kappa = exp(p[2])
  mu = p[3]
  tau_beta = exp(p[4])
  kappa_beta = exp(p[5])
  if(tau>1e-16 && kappa>1e-16 && tau_beta>1e-16 && kappa_beta>1e-16)
  {
    Y_mu = obj$Y - mu
    Q0 = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q.p = kronecker(Diagonal(obj$n.cv-1,1),Q0)
    #Q_beta = tau_beta*Diagonal(777)
    Q_beta  <- tau_beta*(obj$M0.beta + kappa_beta*obj$M1.beta + kappa_beta^2*obj$M2.beta)
    AtQ <- t(obj$X)%*%Q.p
    Rp <- chol(Q_beta  + AtQ%*%obj$X)
    v = solve(t(Rp),AtQ%*%Y_mu)
    l = (obj$n.cv-1)*sum(log(diag(chol(Q0))))-sum(log(diag(Rp))) + sum(log(diag(chol(Q_beta))))  + t(v)%*%v/2 - Y_mu%*%Q.p%*%Y_mu/2
    cat(c(p,-as.double(l[1])),'\n')
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

A = Diagonal(777,1)*quant.Xt[[1]]
for(j in 2:(n.cv-1)){
  A = rBind(A,Diagonal(777,1)*quant.Xt[[j]])
}

obj = list(X=A, A = Diagonal(777*(n.cv-1),1),
           Y = yv, I = Diagonal(777,1), M0 = M0, M1 = M1,M2 = M2,
           N = 777,n.cv = n.cv)
if(use_quant_scaling){
  obj$A = Diagonal(777,1/sqrt(xv))
}

if(smooth.beta){
obj$M0.beta = M0[reo.orig,reo.orig][reo.orig,reo.orig]
obj$M1.beta = M1[reo.orig,reo.orig][reo.orig,reo.orig]
obj$M2.beta = M2[reo.orig,reo.orig][reo.orig,reo.orig]

res <- optim(c(0,0,0,0,0,0), llike.beta,control = list(REPORT = 1,maxit = 20000), obj=obj)
print(paste("res$convergere =",res$convergence))
tau = exp(res$par[1])
kappa = exp(res$par[2])
mu = res$par[3]
tau_beta = exp(res$par[4])
kappa_beta = exp(res$par[5])
Q = tau*kronecker(Diagonal(obj$n.cv-1,1),obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
Q_beta = tau_beta*(obj$M0.beta + kappa_beta*obj$M1.beta + kappa_beta^2*obj$M2.beta)
AtQ <- t(obj$X)%*%Q
Q.post <-  Q_beta + AtQ%*%obj$X

} else {
res <- optim(c(0,0,0,0), llike.mat2,control = list(REPORT = 1,maxit = 20000), obj=obj)
print(paste("res$convergere =",res$convergence))
tau = exp(res$par[1])
kappa = exp(res$par[2])
mu = res$par[3]
sigma2 = exp(res$par[4])
Q = tau*kronecker(Diagonal(obj$n.cv-1,1),
                  obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
Q = obj$A%*%Q%*%obj$A
AtQ <- t(obj$X)%*%Q
Q.post <-  obj$I/sigma2 + AtQ%*%obj$X
}

quant.est <- as.vector(solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu)))
quant.Yspatial = mu + quant.Xe*quant.est


cat(c(mean((exp(quant.Ye) - exp(quant.y))^2),
                 mean((exp(quant.Ye) - exp(quant.Yp))^2),
                 mean((exp(quant.Ye) - exp(quant.Yspatial))^2)))

qplot.yspatial = rep(NA,m*n)
par(mfrow=c(2,2))
qplot.yspatial[obs.ind] = quant.Ye
dim(qplot.yspatial)<- c(m,n)
image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Data")

qplot.yspatial[obs.ind] = quant.Yspatial
dim(qplot.yspatial)<- c(m,n)
image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Prediction")

qplot.yspatial[obs.ind] = quant.est
dim(qplot.yspatial)<- c(m,n)
image.plot(lon.norway,lat.norway,qplot.yspatial, main = "beta")

qplot.yspatial[obs.ind] = quant.Xe
dim(qplot.yspatial)<- c(m,n)
image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Xe")
