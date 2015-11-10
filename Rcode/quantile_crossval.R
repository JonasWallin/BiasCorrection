graphics.off()
rm(list=ls())
load("../data.downscaledERA40.norway.sav")
load("../data/index.norway.1x1kmAggregert.sav")

#############################
## OPTIONS FOR ESTIMATION
#############################

do.plot = TRUE            #Visualise results?
save.plot = FALSE         #Save figure as pdf?
use_log = TRUE           #Use model in log scale?
use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
use_BCM_eval   = TRUE    #Evaluate results on BCM? If false, use ERA40
use_quant_scaling = FALSE #Scale spatial field by covariate?
n.cv = 7                  #Number of cross validation sets
q = 0.95                   #The quantile to do prediction for
smooth.X = FALSE
alpha = 1
season = 4

cat("-- Model details --\n")
cat("Model : ")
if(use_log){ cat("log")}else { cat("linear")}
if(use_quant_scaling){cat(" with ")}else{ cat(" without ")}
cat("quantile scaling \n")
cat("Training Data :")
if(use_BCM_train){cat(" BCM")}else{cat(" ERA40")}
cat("\nEvaluation Data :")
if(use_BCM_eval){cat(" BCM ")}else {cat(" ERA40 ")}
cat("\n")

#############################
## LOAD DATA
#############################
X = data.downscaledERA40.norway
m = dim(X)[1]
n = dim(X)[2]
rm(X)
obs.ind = (index.norway.1x1kmAggregert[,2]-1)*m+index.norway.1x1kmAggregert[,1]

source("climate.utils.R")

source("load.BCM.data_TheNorwayWay.R")

BCM <- cBind(X_train,X_eval)
Y   <- cBind(Y_train, Y_eval)

source("load.climate.data_TheNorwayWay.R")
ERA <- cBind(X_train,X_eval)

ind = 1:2688
if(use_log)
{
  BCM  = log(BCM)
  Y  = log(Y)
  ERA = log(ERA)
}

n <- dim(Y)[2]
ind <- 1:n

n.sub = ceiling(n/n.cv)
ind.sub <- list()

tmp <- matrix(1:(n.sub*n.cv),n.sub,n.cv)
for(i in 1:(n.cv-1)){
  ind.sub[[i]] = tmp[,i]
}
ind.sub[[n.cv]] = intersect(tmp[,n.cv],ind)

quant.BCM <- quant.ERA <- quant.Y <- list()
for(i in 1:n.cv){
  quant.BCM[[i]] <- apply(BCM[,ind.sub[[i]]],1,quantile,probs=c(q))
  quant.ERA[[i]] <- apply(ERA[,ind.sub[[i]]],1,quantile,probs=c(q))
  quant.Y[[i]]   <- apply(Y[,ind.sub[[i]]],1,quantile,probs=c(q))
}

loc = cBind(as.vector(lon.norway),as.vector(lat.norway))[obs.ind,]

oi = rep(FALSE,m*n)
oi[obs.ind] = TRUE
mesh = excursions:::submesh.grid(matrix(oi,m,n),
              list(loc=cBind(as.vector(lon.norway),as.vector(lat.norway))
                          ,dims=c(m,n)))

fem = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2,
                         output = list("c0", "c1", "g1", "g2"))

M0 <- t(fem$g1)%*%Diagonal(mesh$n,1/rowSums(fem$c1))%*%fem$g1
M1 <- fem$g1 + t(fem$g1)
M2 <- fem$c1

idx <- mesh$idx$loc
ridx <- integer(mesh$n)
ridx[idx[!is.na(idx)]] = which(!is.na(idx))
reo <- match(obs.ind,ridx)
ireo = 1:777; ireo[reo] = 1:777
#tmp <- ireo
#ireo <- reo
#reo <- tmp

if(smooth.X){
if(use_log == T){
  load(paste("temp_data/Xsmooth_BCM_log_alpha",alpha,".RData",sep=""))
  quant.BCM <- X_smoothed
  load(paste("temp_data/Xsmooth_ERA_log_alpha",alpha,".RData",sep=""))
  quant.ERA <- X_smoothed
}else{
  load("temp_data/Xsmooth_BCM.RData")
  quant.BCM <- X_smoothed
  load("temp_data/Xsmooth_ERA.RData")
  quant.ERA <- X_smoothed
}
}

for(i in 1:n.cv){
  quant.BCM[[i]] <- quant.BCM[[i]][reo]
  quant.ERA[[i]] <-   quant.ERA[[i]][reo]
  quant.Y[[i]]   <- quant.Y[[i]][reo]
}

if(0){

  llike.mat <- function(p,obj)
  {
    tau = exp(p[1])
    kappa = exp(p[2])
    sigma2 = exp(p[3])
    if(tau>1e-16 && kappa>1e-16 && sigma2>1e-16)
    {
      Q = tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
      Rp = chol(Q + obj$AtA/sigma2)
      v = solve(Rp,obj$AtY/sigma2)
      R = chol(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
      l = log(tau)*obj$n + 2*sum(log(diag(R))) -2*sum(log(diag(Rp))) -obj$N*log(sigma2) + t(v)%*%v -obj$yty/sigma2
      l = as.vector(l)
      return(-l/2)
    } else {
      return(-Inf)
    }
  }

  A = Diagonal(777,1)[reo,]
  A = A[reo,]
  X  = quant.BCM[[i]]
  obj = list(AtA = t(A)%*%A, AtY = as.vector(t(A)%*%X), yty = as.vector(t(X)%*%X),
             M0 = M0, M1 = M1,M2 = M2, n = 777, N = 777)
  cat("smoothing X")
  res <- optim(c(0,0,0), llike.mat,control = list(REPORT = 1), obj=obj)
  tau = exp(res$par[1])
  kappa = exp(res$par[2])
  s2 = exp(res$par[3])

  #kappa = 1
  #s2 = 0.1
  cat(", done \n")
  Q.smoother <- tau*(M0 + kappa*M1 + kappa^2*M2 + t(A)%*%A/s2)
  for(i in 1:n.cv){
    quant.BCM[[i]] <- as.vector(A%*%solve(Q.smoother,t(A)%*%quant.BCM[[i]]/s2))
    quant.ERA[[i]] <- as.vector(A%*%solve(Q.smoother,t(A)%*%quant.ERA[[i]]/s2))
  }
}

qplot = rep(NA,m*n)

if(save.plot){
} else { dev.new() }
par(mfrow = c(3,3))
for(i in 1:n.cv){
  if(save.plot){
    fig_name = paste("../letter/Y",i,".pdf",sep="")
    pdf(fig_name)
  }
  qplot[obs.ind] = quant.Y[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
  if(save.plot){dev.off()}
}

if(save.plot){
} else { dev.new() }

par(mfrow = c(3,3))
for(i in 1:n.cv){
  if(save.plot){
    fig_name = paste("../letter/ERA",i,".pdf",sep="")
    pdf(fig_name)
  }
  qplot[obs.ind] = quant.ERA[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
  if(save.plot){dev.off()}
}

if(save.plot){
} else { dev.new() }

par(mfrow = c(3,3))
for(i in 1:n.cv){
  if(save.plot){
    fig_name = paste("../letter/BCM",i,".pdf",sep="")
    pdf(fig_name)
  }
  qplot[obs.ind] = quant.BCM[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
  if(save.plot){dev.off()}
}


vars <- means <- mse <- matrix(0,n.cv,3)
spatial.coverage <- rep(0,n.cv)
#Run cross validation
for(i in 1:n.cv){
  ind.test <- ind.sub[[i]]
  ind.train <- setdiff(ind,ind.sub[[i]])

  if(use_BCM_train)
  {
    quant.Xt <- quant.BCM[setdiff(1:n.cv,i)]
  }else{
    quant.Xt <- quant.ERA[setdiff(1:n.cv,i)]
  }

  if(use_BCM_eval) {
    quant.Xe <- quant.BCM[[i]]
    name = "BCM"
  } else {
    quant.Xe <- quant.ERA[[i]]
    name = "ERA"
  }
  #Quantiles based on all test data
  quant.y   <- apply(Y[,setdiff(ind,ind.sub[[i]])],1,quantile,probs=c(q))
  quant.y <- quant.y[reo]
  quant.Yt <- quant.Y[setdiff(1:n.cv,i)]
  quant.Ye <- quant.Y[[i]]
  #quant.p <- matrix(0,777,n.cv-1)
  #for(j in 1:(n.cv-1)){
  #  quant.p[,j] <- quant.Yt[[j]]/quant.Xt[[j]]
  #}
  #quant.Yp = quant.Xe*rowMeans(quant.p)
  yv <- quant.Yt[[1]]
  xv <- quant.Xt[[1]]
  xx = cBind(Diagonal(777,1), Diagonal(777,1)*quant.Xt[[1]])
  for(j in 2:(n.cv-1)){
    yv <- c(yv,quant.Yt[[j]])
    xv <- c(xv,quant.Xt[[j]])
    xx = rBind(xx,cBind(Diagonal(777,1),Diagonal(777,1)*quant.Xt[[j]]))
  }

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

  mu_b_ = 0
  llike.mat2 <- function(p,obj)
  {
    tau = exp(p[1])
    kappa = exp(p[2])
    mu = p[3]
    sigma2 = exp(p[4])
    mu_b = mu_b_#p[4]
    #cat(p,"\n")
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
  mu_b = mu_b_#res$par[4]
  Q = tau*kronecker(Diagonal(obj$n.cv-1,1),
                          obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
  Q = obj$A%*%Q%*%obj$A
  AtQ <- t(obj$X)%*%Q
  Q.post <-  obj$I/sigma2 + AtQ%*%obj$X

  quant.est <- as.vector(solve(Q.post,t(obj$X)%*%Q%*%(obj$Y-mu) + mu_b/sigma2))

  quant.Yspatial = mu + quant.Xe*quant.est

  #Compute coverage:
  quant.Yspatial.var <- quant.Xe^2*diag(inla.qinv(Q.post)) +
            diag(inla.qinv(tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)))
  T <- (quant.Ye -quant.Yspatial)/sqrt(quant.Yspatial.var)
  spatial.coverage[i] <- mean(abs(T) < 1.96)

  qplot.yspatial = rep(NA,m*n)
  qplot.yspatial[obs.ind] = quant.Yspatial[ireo]
  dim(qplot.yspatial)<- c(m,n)
  image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Prediction")


  if(use_log){
    vars[i,] <- c(var(exp(quant.Ye) - exp(quant.y)),
                  var(exp(quant.Ye) - exp(quant.Yp)),
                  var(exp(quant.Ye) - exp(quant.Yspatial)))
    means[i,] <- c(mean(exp(quant.Ye) - exp(quant.y)),
                   mean(exp(quant.Ye) - exp(quant.Yp)),
                   mean(exp(quant.Ye) - exp(quant.Yspatial)))
    mse[i,] <- c(mean((exp(quant.Ye) - exp(quant.y))^2),
                 mean((exp(quant.Ye) - exp(quant.Yp))^2),
                 mean((exp(quant.Ye) - exp(quant.Yspatial))^2))

  }else{
    vars[i,] <- c(var(quant.Ye - quant.y),
                  var(quant.Ye - quant.Yp),
                  var(quant.Ye - quant.Yspatial))
    means[i,] <- c(mean(quant.Ye - quant.y),
                   mean(quant.Ye - quant.Yp),
                   mean(quant.Ye - quant.Yspatial))
    mse[i,] <- c(mean((quant.Ye - quant.y)^2),
                 mean((quant.Ye - quant.Yp)^2),
                 mean((quant.Ye - quant.Yspatial)^2))
  }

  cat("Resisual MSE:\n")
  cat("    Historcial data  : ", mse[i,1], "\n")
  cat("    Independent model: ", mse[i,2], "\n")
  cat("    Spatial model    : ", mse[i,3], " (coverage : ",spatial.coverage[i], ")\n")

  if(do.plot){
    qplot.ye = qplot.se = qplot.yspatial = qplot.y = rep(NA,m*n)
    qplot.y[obs.ind] = quant.y[ireo]
    qplot.ye[obs.ind] = quant.Ye[ireo]
    qplot.se[obs.ind] = quant.Yp[ireo]
    qplot.yspatial[obs.ind] = quant.Yspatial[ireo]
    dim(qplot.y)<- dim(qplot.ye)<- dim(qplot.se)<- dim(qplot.yspatial)<- c(m,n)
    dev.new()
    par(mfrow=c(2,2))
    image.plot(lon.norway,lat.norway,qplot.ye, main = "Observation")
    image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Prediction")
    qplot.beta = rep(NA,m*n)
    qplot.beta[obs.ind] = quant.est[ireo]
    dim(qplot.beta) <- c(m,n)
    image.plot(lon.norway,lat.norway,qplot.beta, main = "beta")
    qplot.X = rep(NA,m*n)
    qplot.X[obs.ind] = quant.Xe[ireo]
    dim(qplot.X) <- c(m,n)
    image.plot(lon.norway,lat.norway,qplot.X, main = "X")


    if(save.plot){
      fig_name = paste("../figs/regulaized_prediction_",name,".pdf",sep="")
      pdf(fig_name)
    }else{dev.new()}
    par(mfrow=c(3,3))
    if(use_log == FALSE){
      image.plot(lon.norway,lat.norway,qplot.ye, main = "Observation")
      image.plot(lon.norway,lat.norway,qplot.se,
                  main = paste(name,"Prediction"))
      image.plot(lon.norway,lat.norway,qplot.yspatial,
                  main = paste(name,"Spatial"))
      image.plot(lon.norway,lat.norway,qplot.ye - qplot.y,
                  main ="Resid past data")
      image.plot(lon.norway,lat.norway,qplot.ye - qplot.se,
                  main = "Resid prediction")
      image.plot(lon.norway,lat.norway,qplot.ye - qplot.yspatial,
                  main = "Resid spatial")
      hist(quant.Ye - quant.y,100,
          main= paste("Var = ",formatC(var(quant.Ye - quant.y),digits=2)))
      hist(quant.Ye - quant.Yp,100,
          main= paste("Var = ",formatC(var(quant.Ye - quant.Yp),digits=2)))
      hist(quant.Ye - quant.Yspatial,100,
        main= paste("Var = ",formatC(var(quant.Ye - quant.Yspatial),digits=2)))
    }else{
      image.plot(lon.norway,lat.norway,exp(qplot.ye), main = "Observation")
      image.plot(lon.norway,lat.norway,exp(qplot.se),
                  main = paste(name,"Prediction"))
      image.plot(lon.norway,lat.norway,exp(qplot.yspatial),
                  main = paste(name,"Spatial"))
      image.plot(lon.norway,lat.norway,exp(qplot.ye) - exp(qplot.y),
                  main ="Resid past data")
      image.plot(lon.norway,lat.norway,exp(qplot.ye) - exp(qplot.se),
                  main = "Resid prediction")
      image.plot(lon.norway,lat.norway,exp(qplot.ye) - exp(qplot.yspatial),
                  main = "Resid spatial")
      hist(exp(quant.Ye) - exp(quant.y),100, main= paste("Var = ",
          formatC(var(exp(quant.Ye) - exp(quant.y)),digits=2)))
      hist(exp(quant.Ye) - exp(quant.Yp),100, main= paste("Var = ",
          formatC(var(exp(quant.Ye) - exp(quant.Yp)),digits=2)))
      hist(exp(quant.Ye) - exp(quant.Yspatial),100,main= paste("Var = ",
          formatC(var(exp(quant.Ye) - exp(quant.Yspatial)),digits=2)))
    }
    if(save.plot){dev.off()}
  }
}

