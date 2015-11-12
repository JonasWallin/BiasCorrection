##
# internal function for crossval.R
##
cat('in _crossvalModel12.R')
par_est = solve(t(xx)%*%xx,t(xx)%*%yv)
dim(par_est) <- c(777,2)
alpha_est = par_est[,1]
beta_est = par_est[,2]

quant.Yp0 = alpha_est + beta_est*quant.Xe

xx = Diagonal(777,1)*quant.Xt[[1]]
for(j in 2:(n.cv-1)){
  xx = rBind(xx,Diagonal(777,1)*quant.Xt[[j]])
}
xx = cBind(rep(1,dim(xx)[1]),xx)

par_est = solve(t(xx)%*%xx,t(xx)%*%yv)
alpha_est = par_est[1]
beta_est = par_est[2:778]

quant.Yp = alpha_est + beta_est*quant.Xe

llike.mat2 <- function(p,obj)
{
  tau = exp(p[1])
  kappa = exp(p[2])
  mu = p[3]
  sigma2 = exp(p[4])
  mu_b = 0
  if(tau>1e-16 && kappa>1e-16 && sigma2>1e-16)
  {
    Y_mu = obj$Y - mu
    Q0 = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q.p = kronecker(Diagonal(obj$n.cv-1,1),Q0)
    
    Q.p = obj$A%*%Q.p%*%obj$A
    AtQ <- t(obj$X)%*%Q.p
    Qhat <-  obj$I/sigma2  + AtQ%*%obj$X
    Rp <- chol(Qhat)
    v = solve(t(Rp),AtQ%*%Y_mu+ mu_b/sigma2)
    l = (obj$n.cv-1)*sum(log(diag(chol(Q0))))-sum(log(diag(Rp))) - obj$N*log(sigma2)/2  + t(v)%*%v/2 - Y_mu%*%Q.p%*%Y_mu/2
    l = l - nrow(obj$I)* mu_b^2/(2*sigma2)#correction for mean
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

res <- optim(c(0,0,0,0), llike.mat2,control = list(REPORT = 1,maxit = 20000), obj=obj)
print(paste("res$convergere =",res$convergence))
tau = exp(res$par[1])
kappa = exp(res$par[2])
mu = res$par[3]
sigma2 = exp(res$par[4])
mu_b = 0
Q = tau*kronecker(Diagonal(obj$n.cv-1,1),
                  obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
Q = obj$A%*%Q%*%obj$A
AtQ <- t(obj$X)%*%Q
Q.post <-  obj$I/sigma2 + AtQ%*%obj$X

quant.est <- as.vector(solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu) + mu_b/sigma2))
