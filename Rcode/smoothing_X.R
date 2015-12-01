
if(exists("smooth.X") == FALSE){
  rm(list=ls())
  data_location = '../data/' #NULL
  graphics.off()
  #############################
  ## OPTIONS FOR ESTIMATION
  #############################
  use_BCMs  = c(TRUE,FALSE)
  use_log = TRUE            #Use model in log scale?
  use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
  use_BCM_eval   = TRUE     #Evaluate results on BCM? If false, use ERA40
  use_quant_scaling = FALSE #Scale spatial field by covariate?
  n.cv = 7                  #Number of cross validation sets
  q = 0.95                  #The quantile to do prediction for
  season = 1                #Season to work with (1=winter, 2=spring,...)


  use_BCMs  = c(TRUE, FALSE)
  source('_data_building.R')
}
cat('smoothing_X')


posterior.mean <- function(p, obj)
{
 tau     = exp(p[1])
  kappa   = exp(p[2])
  sigma2 = exp(p[3])
  mu      = p[4]
  Y_mu = obj$Y - mu
  Q0 <- tau*(obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 )
  Q <- kronecker(Diagonal(obj$n.cv,1),Q0)
  Q.post <- Q + t(obj$A)%*%obj$A/sigma2
  Rp <- chol(Q.post)
  b <- t(obj$A)%*%Y_mu/sigma2
  v = solve(t(Rp),b)

  xv  <- as.vector(obj$A%*%(mu+solve(Rp, v)))
  dim(xv) <- c(777,obj$n.cv)
  X_smooth <- list()
  for(i in 1:obj$n.cv){
    X_smooth[[i]] = xv[,i]
  }

  return(X_smooth)
}

llike.mat2 <- function(p,obj)
{
  tau     = exp(p[1])
  kappa   = exp(p[2])
  sigma2 = exp(p[3])
  mu      = p[4]
  if(tau>1e-16 && kappa>1e-16)
  {
    Y_mu = obj$Y - mu
    Q0 <- tau*(obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 )
    Q <- kronecker(Diagonal(obj$n.cv,1),Q0)
    Q.post <- Q + obj$AtA/sigma2
    Rp <- chol(Q.post)
    b <- t(obj$A)%*%Y_mu/sigma2
    v = solve(t(Rp),b)
    l = obj$n.cv*sum(log(diag(chol(Q0)))) -sum(log(diag(Rp))) - length(Y_mu)*log(sigma2)/2 + t(v)%*%v/2 - t(Y_mu)%*%Y_mu/(2*sigma2)
    #cat(p,-as.double(l[1]),'\n')
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

llike.cov <- function(p,obj)
{
  cat(p)

  n <- length(obj$Y[[1]])
  n.rep <- length(obj$Y)
  mu      = p[5]
  Sigma.X <- matern.cov(obj$D,kappa=exp(p[2]),sigma2 = exp(p[1]),nu=exp(p[3]))
  Sigma.E <- exp(p[4])*Diagonal(n)

  Sigma <- Sigma.X + Sigma.E
  R <- try(chol(Sigma),silent=TRUE)
  if(class(R) == "try-error")
    return(-Inf)

  l <- -n.rep*sum(log(diag(R)))
  for(i in 1:n.rep){
    v <- solve(t(R),obj$Y[[i]]-mu)
    l <- l - t(v)%*%v/2
  }
  cat(' ',-as.double(l[1]),'\n')
  return(-as.double(l[1]))
}


posterior.mean.cov <- function(p, obj)
{
  n <- length(obj$Y[[1]])
  n.rep <- length(obj$Y)
  mu      = p[5]
  Sigma.X <- matern.cov(obj$D,kappa=exp(p[2]),sigma2 = exp(p[1]),nu=exp(p[3]))
  sigma2 = exp(p[4])
  Sigma.post = solve(solve(Sigma.X) + Diagonal(n)/sigma2)
  X_smooth <- list()
  for(i in 1:n.rep){
    X_smooth[[i]] = mu + Sigma.post%*%(obj$Y[[i]]-mu)/sigma2
  }
  return(X_smooth)
}

for(kk in 1:length(use_BCMs)){

  use_BCM <- use_BCMs[kk]
  if(use_BCM)
  {
    quant.Xt <- quant.BCM
  }else{
    quant.Xt <- quant.ERA
  }

  n.cv = length(quant.Xt)
  if(use.cov){

    obj <- list(D=as.matrix(dist(loc)),Y = quant.Xt)
    if(use_BCM){
      p0 <- c(-0.851345, -2.389743, -0.7534426, -2.207971, 2.217278)
    } else {
      p0 <- c(-1.268291, -1.046623, -0.3034919, -2.552335, 2.359651)
    }
    res <- optim(p0, llike.cov,
                 control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
    X_smooth <- posterior.mean.cov(res$par,obj)
  } else {
    A = kronecker(Diagonal(n.cv),Ax)#Diagonal(777*n.cv,1)
    xx = quant.Xt[[1]]
    for(j in 2:n.cv){
      xx = cBind(xx,quant.Xt[[j]])
    }
    dim(xx) <- c(777*n.cv,1)

    obj <- list(M0 = M0,#[reo.orig,reo.orig][reo.orig,reo.orig],
              M1 = M1,#[reo.orig,reo.orig][reo.orig,reo.orig],
              M2 = M2,#[reo.orig,reo.orig][reo.orig,reo.orig],
              A = A,
              AtA = t(A)%*%A,
              n.cv = length(quant.Xt),
              Y = xx)

    p0 <- c(0.1817722, -1.172809, -2.527773, 2.351842)
    res <- optim(p0, llike.mat2,
                  control = list(REPORT = 1,maxit = 20000), obj=obj)
    print(paste("res$convergere =",res$convergence))
    X_smooth <- posterior.mean(res$par, obj)
  }

  if(exists("smooth.X") == FALSE){
    file_name <- paste(data_location,"Xsmooth_",sep="")

    if(use_BCM)
    {
      file_name <- paste(file_name,"BCM",sep="")
    }else{
      file_name <- paste(file_name,"ERA",sep="")
    }

    if(use_log)
      file_name <- paste(file_name,"_log",sep="")

    file_name <- paste(file_name,"_season_", season, sep = "")
    file_name <- paste(file_name,".RData", sep="")

    save(X_smooth,file=file_name)
  }
}


cat('done smoothing\n')
