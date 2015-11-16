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
  if(tau>1e-16 && kappa>1e-16 && tau_beta>1e-16 && kappa_beta>1e-16)
  {
    Y_mu = obj$Y - mu
    Q0 = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    Q.p = kronecker(Diagonal(obj$n.cv-1,1),Q0)
    #Q_beta = tau_beta*Diagonal(777)
    Q_beta  <- tau_beta*(obj$M0 + kappa_beta*obj$M1 + kappa_beta^2*obj$M2)
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
#Ai = inla.spde.make.A(mesh,loc)[reo,] #Ai*v = v[ireo]
Ai = Diagonal(777)[reo,][reo,]
#Ai  =Ai%*%Ai
A = A%*%Ai

obj = list(X=A, Y = yv, M0 = M0, M1 = M1,M2 = M2, N = 777,n.cv = n.cv)

p0 = c(-3.37165, 4.06397, 2.504563, 1.054453, -0.2670174)
res <- optim(p0, llike.mat2,control = list(REPORT = 1,maxit = 20000), obj=obj)
print(paste("res$convergere =",res$convergence))
tau = exp(res$par[1])
kappa = exp(res$par[2])
mu = res$par[3]
tau_beta = exp(res$par[4])
kappa_beta = exp(res$par[5])
Q = tau*kronecker(Diagonal(obj$n.cv-1,1),obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
Q_beta = tau_beta*(obj$M0 + kappa_beta*obj$M1 + kappa_beta^2*obj$M2)
AtQ <- t(obj$X)%*%Q
Q.post <-  Q_beta + AtQ%*%obj$X

beta <- as.vector(Ai%*%solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu)))
quant.Yspatial = as.vector(mu + quant.Xe*beta)

cat(c(mean((exp(quant.Ye) - exp(quant.y))^2),
                 mean((exp(quant.Ye) - exp(quant.Yp))^2),
                 mean((exp(quant.Ye) - exp(quant.Yspatial))^2)))

qplot.yspatial = rep(NA,m*n)
qplot.yspatial[obs.ind] = quant.Xe
dim(qplot.yspatial)<- c(m,n)
image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Prediction")
