##
# internal function for crossval.R
##

llike.cov <- function(p,obj)
{
  n <- length(obj$Y[[1]])
  n.rep <- length(obj$Y)

  if(smooth.beta){
    n.beta <- 2
    Sigma.beta <- matern.cov(obj$D,kappa=exp(p[2]),
                             sigma2 = exp(p[1]),nu=obj$nu)
  } else {
    n.beta <- 1
    Sigma.beta <- exp(p[1])*Diagonal(n)
  }
  if(smooth.error){
    n.error <- 2
    Sigma.E <- matern.cov(obj$D,kappa=exp(p[n.beta+2]),
                          sigma2 = exp(p[n.beta+1]),nu=obj$nu)
  } else {
    n.error <- 1
    Sigma.E <- exp(p[n.beta+1])*Diagonal(n)
  }
  mu = p[n.beta+n.error+1]

  l <- 0
  for(ii in 1:n.rep){
    X <- Diagonal(777,obj$X[[ii]])
    Sigma <- X%*%Sigma.beta%*%t(X) + Sigma.E
    Sigma = (Sigma + t(Sigma))/2
    R <- chol(Sigma)
    v <- solve(R,obj$Y[[ii]]-mu)
    l <- l - sum(log(diag(R))) - t(v)%*%v/2
  }
  cat(c(p,-as.double(l[1])),'\n')
  return(-as.double(l[1]))
}


llike.Esmooth <- function(p,obj)
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

    AtQ <- t(obj$X)%*%Q.p
    Qhat <-  obj$I/sigma2  + AtQ%*%obj$X
    Rp <- chol(Qhat)
    v = solve(t(Rp),AtQ%*%Y_mu)
    l = (obj$n.cv-1)*sum(log(diag(chol(Q0))))-sum(log(diag(Rp))) - dim(obj$I)[1]*log(sigma2)/2  + t(v)%*%v/2 - Y_mu%*%Q.p%*%Y_mu/2

    #l = l + dnorm(sigma2,var(yv)/var(xv),0.01,log=TRUE)
    #cat(c(p,-l[1]),'\n')
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

llike.Bsmooth <- function(p,obj)
{
  tau     = exp(p[1])
  kappa   = exp(p[2])
  sigma2 = exp(p[3])
  mu      = p[4]
  if(tau>1e-16 && kappa>1e-16)
  {
    Y_mu = obj$Y - mu
    Q <- tau*(obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 )
    Q.post <- Q + obj$XtX/sigma2
    Rp <- chol(Q.post)
    b <- t(obj$X)%*%Y_mu/sigma2
    v = solve(t(Rp),b)
    l = sum(log(diag(chol(Q)))) -sum(log(diag(Rp))) - length(Y_mu)*log(sigma2)/2 + t(v)%*%v/2 - t(Y_mu)%*%Y_mu/(2*sigma2)
    #cat(p,-as.double(l[1]),'\n')
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

llike.BsmoothEsmooth <- function(p,obj)
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
    #l = l + dnorm(p[2],20,0.1,log=TRUE)
    cat(c(p,-as.double(l[1])),'\n')
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

matern.cov <- function(D,kappa,sigma2,nu)
{
  Sigma <- sigma2*2^(1-nu)*(kappa*D)^nu * besselK(kappa*D,nu)/gamma(nu)
  if(is.matrix(Sigma)){
    diag(Sigma) <- sigma2
  } else {
    Sigma[1] <- sigma2
  }
  return(Sigma)
}

if(use.cov) {

  obj = list(X=quant.Xt,Y = quant.Yt,D=as.matrix(dist(loc)),nu=0.5)
  if(smooth.beta){ n.beta = 2}else {n.beta = 1}
  if(smooth.error) {n.error = 2}else {n.error = 1}
  p0 = rep(0,n.beta+n.error+1)
  if(smooth.beta ==TRUE &&  smooth.error == FALSE){
    p0 = c(-3.540655, -0.6742775, -4.902866, 2.302089)
  }
  res <- optim(p0, llike.cov,
               control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
  n <- length(obj$Y[[1]])
  n.rep <- length(obj$Y)

  if(smooth.beta){
    Sigma.beta <- matern.cov(obj$D,kappa=exp(res$par[2]),
                             sigma2 = exp(res$par[1]),nu=obj$nu)
  } else {
    n.beta <- 1
    Sigma.beta <- exp(res$par[1])*Diagonal(n)
  }
  if(smooth.error){
    n.error <- 2
    Sigma.E <- matern.cov(obj$D,kappa=exp(res$par[n.beta+2]),
                          sigma2 = exp(res$par[n.beta+1]),nu=obj$nu)
  } else {
    n.error <- 1
    Sigma.E <- exp(res$par[n.beta+1])*Diagonal(n)
  }
  mu = res$par[n.beta+n.error+1]
  Q.E = solve(Sigma.E)
  Q.post = solve(Sigma.beta)
  b = 0
  for(ii in 1:n.rep){
    X <- Diagonal(777,obj$X[[ii]])
    Q.post = Q.post + X%*%Q.E%*%X
    b = b + X%*%Q.E%*%(obj$Y[[ii]]-mu)
  }
  quant.est = as.double(solve(Q.post,b))

} else {
  A = Diagonal(777,1)*quant.Xt[[1]]
  for(j in 2:(n.cv-1)){
    A = rBind(A,Diagonal(777,1)*quant.Xt[[j]])
  }

  obj = list(X=A,Y = yv, I = Diagonal(777,1),n.cv = n.cv,
           XtX = t(A)%*%A,
           M0 = M0[ireo.orig,ireo.orig][ireo.orig,ireo.orig],
           M1 = M1[ireo.orig,ireo.orig][ireo.orig,ireo.orig],
           M2 = M2[ireo.orig,ireo.orig][ireo.orig,ireo.orig])

  if(smooth.beta && smooth.error)
  {
    obj$M0.beta = M0
    obj$M1.beta = M1
    obj$M2.beta = M2

    res <- optim(c(0,0,0,0,0), llike.beta,
               control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
    tau = exp(res$par[1])
    kappa = exp(res$par[2])
    mu = res$par[3]
    tau_beta = exp(res$par[4])
    kappa_beta = exp(res$par[5])
    Q = tau*kronecker(Diagonal(obj$n.cv-1,1),
                      obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q_beta =tau_beta*(obj$M0.beta+kappa_beta*obj$M1.beta+
                    kappa_beta^2*obj$M2.beta)
    AtQ <- t(obj$X)%*%Q
    Q.post <-  Q_beta + AtQ%*%obj$X
    quant.est <- as.vector(solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu)))
  } else if (smooth.error) {
    res <- optim(c(0,0,0,0), llike.Esmooth,
               control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
    tau = exp(res$par[1])
    kappa = exp(res$par[2])
    mu = res$par[3]
    sigma2 = exp(res$par[4])
    Q = tau*kronecker(Diagonal(obj$n.cv-1,1),
                  obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    AtQ <- t(obj$X)%*%Q
    Q.post <-  obj$I/sigma2 + AtQ%*%obj$X
    quant.est <- as.vector(solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu)))
  } else if (smooth.beta) {
    res <- optim(c(0,0,0,0), llike.Bsmooth,
                control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
    tau = exp(res$par[1])
    kappa = exp(res$par[2])
    mu = res$par[3]
    sigma2 = exp(res$par[4])

    Q = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q.post <-  Q + obj$XtX/sigma2
    quant.est <- as.vector(solve(Q.post,t(obj$X)%*%(obj$Y-mu)/sigma2))
  }
}

quant.Yspatial = as.double(mu + quant.Xe*quant.est)


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
