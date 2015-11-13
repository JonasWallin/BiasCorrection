
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
  tau_eps = exp(p[3])
  mu      = p[4]
  
  Q_X  <- tau    * (obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 )
  
  n <- dim(obj$M0)[1]
  N <- length(obj$X_loc)
  QA_xt    <- tau_eps * t(obj$A_x)
  Q_hat <- Q_X + QA_xt%*%obj$A_x 
  Rp <- chol(Q_hat)
  b_0 <- Q_X%*%rep(mu, n)
  yty <- 0
  vtv <- 0
  X_smooth <- obj$X_loc
  for(i in 1:N){
    Y <- obj$X_loc[[i]]
    b    <- QA_xt%*%Y + b_0
    v    <- solve(t(Rp), b)
    X_smooth[[i]]  <- as.vector(solve(Rp, v))
    
  }
  return(X_smooth)
}

llike.mat2 <- function(p,obj)
{
  tau     = exp(p[1])
  kappa   = exp(p[2])
  tau_eps = exp(p[3])
  mu      = p[4]
  if(tau>1e-16 && kappa>1e-16 && tau_eps>1e-16 )
  {
    
    
    Q_X  <- tau    * (obj$M0 + kappa * obj$M1 + kappa^2*obj$M2 )
    
    n <- dim(obj$M0)[1]
    N <- length(obj$X_loc)
    QA_xt    <- tau_eps * t(obj$A_x)
    Q_hat <- Q_X + QA_xt%*%obj$A_x 
    Rp <- chol(Q_hat)
    b_0 <- Q_X%*%rep(mu, n)
    yty <- 0
    vtv <- 0
    for(i in 1:N){
      Y <- obj$X_loc[[i]]
      b    <- QA_xt%*%Y + b_0
      v    <- solve(t(Rp), b)
      #v   <- solve(Rp, v)
      yty <- yty +  tau_eps * (t(Y)%*%Y)/2 
      vtv <- vtv +  (t(v)%*%v)/2
    }
    
    #determinants
    l = - N * sum(log(diag(Rp))) + N * sum(log(diag(chol(Q_X))))  + 0.5 * N * n * log(tau_eps)
    #
    l = l + vtv - yty
    l = l - N * (t(rep(mu, n))%*% b_0)/2#correction for mean
    return(-as.double(l[1]))
  } else {
    return(-Inf)
  }
}

for(kk in 1:length(use_BCMs)){
    
  use_BCM <- use_BCMs[kk]
  if(use_BCM)
  {
    quant.Xt <- quant.BCM
  }else{
    quant.Xt <- quant.ERA
  }
  
   obj <- list(M0 = M0, M1 = M1, M2 = M2, reo = reo, ireo = ireo)
   obj$A_x         <-  Diagonal(length(quant.Xt[[1]]))
   obj$X_loc       <- quant.Xt 


  res <- optim(c(0,0,0,0), llike.mat2,control = list(REPORT = 1,maxit = 20000), obj=obj)
  print(paste("res$convergere =",res$convergence))
  X_smooth <- posterior.mean(res$par, obj)
  
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
