###
# This program runs MCMC for beta, and parameters using a previously smoothed X as a fixed covariates
# If running this file make sure that the parameters (like n.cv) matches the smoothing files loaded
# build new files using smoothning_X
###

graphics.off()
rm(list=ls())
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
FILE_PATH <- dirname(frame_files[[length(frame_files)]])
setwd('../')
load("../data.downscaledERA40.norway.sav")
load("../data/index.norway.1x1kmAggregert.sav")

#############################
## OPTIONS FOR ESTIMATION
#############################

kappa_beta_fixed <- FALSE


do.plot = TRUE            #Visualise results?
save.plot = FALSE         #Save figure as pdf?
save.res  <- TRUE        # store data
use_log = TRUE           #Use model in log scale?
use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
use_BCM_eval   = TRUE    #Evaluate results on BCM? If false, use ERA40
use_quant_scaling = FALSE #Scale spatial field by covariate?
n.cv = 7                  #Number of cross validation sets
q = 0.95                   #The quantile to do prediction for
sim      = 10000
burning <- 5000
alpha      <- 1
kappa_y    <- 1
kappa_beta <- 1
tau_y      <- 1
tau_beta    <- 0
mu_beta     <- 0
mu_y        <- 0
season      <- 3

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
source("mcmc_sampling.R")
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
#  quant.BCM[[i]] <- apply(BCM[,ind.sub[[i]]],1,quantile,probs=c(q))
#  quant.ERA[[i]] <- apply(ERA[,ind.sub[[i]]],1,quantile,probs=c(q))
  quant.Y[[i]]   <- apply(Y[,ind.sub[[i]]],1,quantile,probs=c(q))
}



if(0){
  qplot = rep(NA,m*n)
  x11()
  par(mfrow = c(3,3))
  for(i in 1:n.cv){
    qplot[obs.ind] = quant.Y[[i]]; dim(qplot) <- c(m,n)
    image.plot(lon.norway,lat.norway,qplot)
  }

  x11()
  par(mfrow = c(3,3))
  for(i in 1:n.cv){
    qplot[obs.ind] = quant.ERA[[i]]; dim(qplot) <- c(m,n)
    image.plot(lon.norway,lat.norway,qplot)
  }

  x11()
  par(mfrow = c(3,3))
  for(i in 1:n.cv){
    qplot[obs.ind] = quant.BCM[[i]]; dim(qplot) <- c(m,n)
    image.plot(lon.norway,lat.norway,qplot)
  }
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

coverage <- rep(0,n.cv)

for(i in 1:n.cv){
 # quant.BCM[[i]] <- quant.BCM[[i]]   #[reo]
#  quant.ERA[[i]] <-   quant.ERA[[i]] #[reo]
  quant.Y[[i]]   <- quant.Y[[i]]     #[reo]
}

quant.Yt <- quant.Y[setdiff(1:n.cv,i)]
vars <- matrix(0,n.cv,3)
#Run cross validation

quant_samples <- list()
Results <- matrix(0, nrow= n.cv, ncol = 4)
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
  quant.Ye <- quant.Y[[i]]
  quant.Yspatial2 <- rep(0,length(quant.Ye))

  quant_samples[[i]] <- matrix(0, nrow= sim - burning +1, ncol = length(quant.Ye) )
 if(i == 1){
   quant.Yspatial2_vec <- matrix(0,nrow = length(quant.Ye), ncol = n.cv)

  obj <- list(M0 = M0, M1 = M1, M2 = M2, reo = reo, ireo = ireo)
  obj$kappa_eps   <- log(kappa_y)
  obj$tau_eps     <- tau_y
  obj$tau_beta    <- (tau_beta)
  obj$kappa_beta  <- log(kappa_beta)
  obj$mu_beta     <- mu_beta
  obj$mu_y        <- mu_y
  #priors beta#
  obj$tau_beta_beta  <- 10^-16
  obj$tau_beta_alpha <- 1
  obj$kappa_beta_lambda <- 1
  #priors eps#
  obj$tau_eps_beta  <- 10^-16
  obj$tau_eps_alpha <- 1
  obj$kappa_eps_lambda <- 1
  #
  ###
  obj$A_x         <-  Diagonal(length(quant.Y[[1]]))[,reo]#Diagonal(length(quant.Y[[1]])) #[reo,]
  #obj$A_x         <-  obj$A_x[reo,]
  obj$A_eps       <- Diagonal(length(quant.Y[[1]]))[,reo]
  obj$A_beta      <- Diagonal(length(quant.Y[[1]]))[,reo]
  obj$beta        <- rep(1,length(quant.Xt[[1]]))


 }
 obj$X           <- quant.Xt
  for(j in 1:length(obj$X))
   {  obj$X[[j]]           <- quant.Xt[[j]][ireo] }

 obj$Y           <- quant.Yt
 #adaptive parameters
 obj$MH_beta <- setup_AMCMC_RR()
 obj$MH_y    <- setup_AMCMC_RR()



  kappa_eps_vec <- rep(0,sim)
  tau_eps_vec   <- rep(0,sim)
  mu_eps_vec    <- rep(0,sim)
kappa_beta_vec <- rep(0,sim)
  tau_beta_vec   <- rep(0,sim)
  mu_beta_vec    <- rep(0,sim)
  beta_vec_mean     <- rep(0,length(obj$beta))
  kappa_x_mean      <- 0
  tau_x_eps_mean    <- 0
  tau_x_mean        <- 0
  mu_x_mean         <- 0
  mu_eps_mean       <- 0
 cat("\n")
  for(j in 1:sim){
    if(j %% 100 == 1){
    cat("*",sep="",append=T)
    }
    if(j %% 1000 == 0){
      cat("+",sep="",append=T)
    }
      #cat('j = ',j,'\n')
      # cat('sampling:  ')



     # cat('Y, ')
      # sampling  \Theta_Y |Y, \beta, X
      res       <- sample_theta_eps(obj)
      kappa_eps_vec[j] <- res$kappa
      tau_eps_vec[j]   <- res$tau
      mu_eps_vec[j]    <- res$mu
      obj$kappa_eps  <- res$kappa
      obj$tau_eps    <- res$tau
      obj$mu_y       <- res$mu
      obj$MH_y       <- AMCMC_RR(obj$MH_y, res$accept)


     # cat('beta\n')
      # sampling \beta \Theta_\beta | Y, \beta, X

      obj$beta  <- sample_beta(obj)

      if(kappa_beta_fixed)
      {
        res       <- sample_theta_beta_kappa_fixed(obj)
      }else{
        res       <- sample_theta_beta(obj)
      }
      kappa_beta_vec[j] <- res$kappa
      tau_beta_vec[j]   <- res$tau
      mu_beta_vec[j]    <- res$mu
      obj$kappa_beta <- res$kappa
      obj$tau_beta   <- res$tau
      obj$mu_beta    <- res$mu
      if(kappa_beta_fixed==FALSE)
        obj$MH_beta <- AMCMC_RR(obj$MH_beta, res$accept)

      if(j >= burning)
      {
        beta_vec_mean     <- obj$beta + beta_vec_mean
        mu_eps_mean       <- mu_eps_mean  + obj$mu_y
        mu.post <- obj$mu_y + obj$A_x%*%quant.Xe[ireo]*(obj$A_beta%*%obj$beta)
        Q.post <- obj$tau_eps*(obj$M0 + exp(obj$kappa_eps)*obj$M1 + exp(2* obj$kappa_eps)*obj$M2)
        samp <- inla.qsample(n=1,mu=mu.post,Q=Q.post)
        if(use_log){
          quant.Yspatial2 <- quant.Yspatial2+ as.vector(exp(samp))
          quant_samples[[i]][j - burning +1, ] <- as.vector(exp(samp))
        }else{
          quant.Yspatial2 <- quant.Yspatial2+ as.vector(samp)
          quant_samples[[i]][j - burning +1, ] <- as.vector(samp)
        }
        #qplot = rep(NA,m*n)
        #qplot[obs.ind] =  as.vector((quant_samples[[i]][j - burning +1,]))
        #dim(qplot)<- c(m,n)
        #x11()
        #image.plot(lon.norway,lat.norway,qplot, main = "Prediction")
      }
  }
    beta_vec_mean   = beta_vec_mean / (sim - burning + 1)
    mu_eps_mean     = mu_eps_mean /(sim - burning + 1)
    quant.Yspatial =  as.vector(mu_eps_mean + (obj$A_x%*%quant.Xe[ireo]) * (obj$A_beta%*%beta_vec_mean))
    Y_mcmc_quantiles <- apply(quant_samples[[i]],2, function(x){ quantile(x,c(0.025,0.5,0.975))})
    quant.Yspatial2 <-  Y_mcmc_quantiles[2,]
    quant.Yspatial2_vec[,i] <- quant.Yspatial2

  if(0){
  x11()
  par(mfrow=c(1, 3))
  qplot.yspatial = rep(NA,m*n)
  qplot.yspatial[obs.ind] =  as.vector((quant.Xe))
  dim(qplot.yspatial)<- c(m,n)
  image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Prediction")


  qplot.yspatial[obs.ind] = quant.Ye
  image.plot(lon.norway,lat.norway,qplot.yspatial, main = "Y_true")

  qplot.yspatial[obs.ind] = as.vector(quant.Yspatial) - quant.Ye
  image.plot(lon.norway,lat.norway,qplot.yspatial, main = "residuals")
  }
  if(use_log){
    var_pred  <-  var( exp(quant.Ye) - quant.Yspatial2)
    mean_pred <- mean( exp(quant.Ye) - quant.Yspatial2)
  }else{
    var_pred  <-  var( quant.Ye - quant.Yspatial2)
    mean_pred <- mean( quant.Ye - quant.Yspatial2)
  }
  MSE <- var_pred + mean_pred^2
  
  cat("\n v2: var_",i," = ", var_pred," mean = ", mean_pred,sep="")
  coverage[i] <- mean((exp(quant.Ye)< Y_mcmc_quantiles[3,])*(exp(quant.Ye)>= Y_mcmc_quantiles[1,]))
  cat("\n     coverage_",i," = ",coverage[i], sep='')
    
 
    
  Results[i, ] <- c(mean_pred, var_pred, MSE, coverage[i])
}

x11()
par(mfrow=c(1, 3))
plot(kappa_eps_vec[10:sim])
plot(tau_eps_vec[10:sim])
plot(mu_eps_vec[10:sim])

x11()
par(mfrow=c(1, 3))
plot(kappa_beta_vec[10:sim])
plot(tau_beta_vec[10:sim])
plot(mu_beta_vec[10:sim])

qplot = rep(NA,m*n)
if(save.plot){
} else { x11() }
par(mfrow = c(3,3))
for(i in 1:n.cv){
  if(save.plot){
    fig_name = paste("../letter/pred_",i,".pdf",sep="")
    pdf(fig_name)
  }
  qplot[obs.ind] = quant.Yspatial2_vec[,i]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
  if(save.plot){dev.off()}


}

if(save.plot)
{

  use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
  use_BCM_eval   = TRUE    #Evaluate results on BCM? If false, use ERA40
  filename = "predict"
  if(use_BCM_eval)
  {
    filename = paste(filename,"_BCM",sep="")
  }else{
    filename = paste(filename,"_ERA",sep="")

  }
  save(quant.Yspatial2_vec, file = paste(filename,".RData",sep=""))
}
setwd(FILE_PATH)

if(save.res){
  
  filename = "Result/crossval_"
  if(use_BCM_eval){
    filename = paste(filename,"_BCM",sep="")
  }else{
    filename = paste(filename,"_ERA",sep="")
  }
  filename = paste(filename,'_',season,".txt",sep="")
  dimnames(Results)[[2]] <-list('mean','var','mse','coverage')
  write.table(Results, file = filename,row.names = F)
}

