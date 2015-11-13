##
# internal function for crossvalModel3.R
##




llike.mat2 <- function(p,obj)
{
  tau = exp(p[1])
  kappa = exp(p[2])
  mu = p[3]
  tau_beta = exp(p[4])
  kappa_beta = exp(p[5])
  mu_b = rep(0, obj$N)
  if(tau>1e-16 && kappa>1e-16 && tau_beta>1e-16 && kappa_beta > 1e-16)
  {
    
    
    Q_eps   <- obj$A_eps%*%(tau    * (obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 ))%*%t(obj$A_eps)
    Q_beta  <- tau_beta   * (obj$M0  + kappa_beta   * obj$M1 + kappa_beta^2*obj$M2)
    
    n <- dim(obj$M0)[1]
    
    
    Q_hat <- Q_beta
    b <- Q_beta%*%rep(obj$mu_beta ,n)
    yQy <- 0
    for( i in 1:length(obj$X)){
      A_x_i = Diagonal(length(obj$A_x%*%obj$X[[1]]),x = as.vector(obj$A_x%*%obj$X[[i]]))%*%obj$A_beta
      QA_x_it <- t(A_x_i)%*%Q_eps
      Q_hat <- Q_hat  + QA_x_it%*%A_x_i
      b    <- b + QA_x_it%*%(obj$Y[[i]] - mu)
      Y_mu <- (obj$Y[[i]] - mu)
      yQy <- yQy +  t(Y_mu)%*%Q_eps%*%Y_mu/2 
      
    }

    
    Rp <- chol(Q_hat)
    v = solve(t(Rp), b)
    l = (obj$n.cv-1)*sum(log(diag(chol(Q_eps)))) + sum(log(diag(chol(Q_beta))))  -sum(log(diag(Rp))) 
    l = l + t(v)%*%v/2 - yQy
    #l = l - (mu_b%*% Q_beta %*%mu_b)/2#correction for mean
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}



obj <- list(M0 = M0, M1 = M1, M2 = M2, reo = reo, ireo = ireo, n.cv = n.cv, N = 777, mu_beta = 0)
obj$A_x         <-  Diagonal(length(quant.Y[[1]]))
#obj$A_x         <-  obj$A_x[reo,]
obj$A_eps       <- Diagonal(length(quant.Y[[1]]))
obj$A_beta      <- Diagonal(length(quant.Y[[1]]))
obj$X           <- quant.Xt
for(j in 1:length(obj$X))
{  obj$X[[j]]           <- quant.Xt[[j]] }

obj$Y           <- quant.Yt

res <- optim(c(0,0,0,0,0), llike.mat2,control = list(REPORT = 1,maxit = 20000), obj=obj)
print(paste("res$convergere =",res$convergence))
tau = exp(res$par[1])
kappa = exp(res$par[2])
mu = res$par[3]
tau_beta = exp(res$par[4])
kappa_beta = exp(res$par[5])
mu_b = rep(0, obj$N)
Q_eps   <- obj$A_eps%*%(tau    * (obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 ))%*%t(obj$A_eps)
Q_beta  <- tau_beta   * (obj$M0  + kappa_beta   * obj$M1 + kappa_beta^2*obj$M2)

n_ <- dim(obj$M0)[1]


Q.post <- Q_beta
b <- Q_beta%*%rep(obj$mu_beta ,n_)
for( j in 1:length(obj$X)){
  A_x_i = Diagonal(length(obj$A_x%*%obj$X[[1]]),x = as.vector(obj$A_x%*%obj$X[[j]]))%*%obj$A_beta
  QA_x_it <- t(A_x_i)%*%Q_eps
  Q.post <- Q.post  + QA_x_it%*%A_x_i
  b    <- b + QA_x_it%*%(obj$Y[[j]] - mu)
  
}
Rp   <- chol(Q.post)
v    <- solve(t(Rp),b)
beta <- solve(Rp, v)
A_beta = Diagonal(length(obj$A_beta%*%beta),as.vector(obj$A_beta%*%beta))%*%obj$A_x
