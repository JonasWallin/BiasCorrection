graphics.off()
rm(list=ls())
###
# Data location, if one has downloaded data the files needed
# spesifiy data location below (if NULL will download the files)
###
data_location = '../data/' #NULL


require(excursions)
library(Matrix)
library(fields)
library(INLA)

#############################
## OPTIONS FOR ESTIMATION
#############################

do.plot = TRUE            #Visualise results?
use_log = TRUE            #Use model in log scale?
use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
use_BCM_eval   = TRUE     #Evaluate results on BCM? If false, use ERA40
n.cv = 7                  #Number of cross validation sets
q = 0.95                  #The quantile to do prediction for
smooth.X = TRUE           #Smooth the covariate?
smooth.beta = TRUE       #Smooth prior for beta?
smooth.error = TRUE       #Smooth prior for eps?
season = 1                #Season to work with (1=winter, 2=spring,...)
alpha = 2                 #Smoothness of random fields
use.cov = TRUE
source('_data_building.R')

if(smooth.X){
  quant.BCM <- load_smooth('BCM', data_location)
  quant.ERA <- load_smooth('ERA', data_location)
}

for(i in 1:n.cv){
  quant.ERA[[i]] = as.double(quant.ERA[[i]])
  quant.BCM[[i]] = as.double(quant.BCM[[i]])
}

if(0){
file_name = paste(data_location,"ERA")
if(use_log){ file_name <- paste(file_name,"_log",sep="")}
file_name <- paste(file_name,"_season_", season,".RData", sep="")
X_smooth <- quant.ERA
save(X_smooth,file=file_name)
file_name = paste(data_location,"BCM")
if(use_log){ file_name <- paste(file_name,"_log",sep="")}
file_name <- paste(file_name,"_season_", season,".RData", sep="")
X_smooth <- quant.BCM
save(X_smooth,file=file_name)
}
###########################
## Plot data
###########################
#ggplot() +  geom_point(aes(x=loc[,1],y=loc[,2],colour=quant.Xt[[1]]), size=2,  alpha=1) + scale_colour_gradientn(colours=tim.colors(100))

if(0){
m = 58
n = 63

qplot = rep(NA,m*n)
dev.new()
par(mfrow = c(3,3))
for(i in 1:n.cv){
  qplot[obs.ind] = quant.Y[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
}

dev.new()
par(mfrow = c(3,3))
for(i in 1:n.cv){
  qplot[obs.ind] = quant.ERA[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
}


dev.new()
par(mfrow = c(3,3))
for(i in 1:n.cv){
  qplot[obs.ind] = quant.BCM[[i]][ireo]; dim(qplot) <- c(m,n)
  image.plot(lon.norway,lat.norway,qplot,ylab="",xlab="")
}
}
##########################
## Run cross validation
##########################
vars <- means <- mse <- matrix(0,n.cv,3)
spatial.coverage <- rep(0,n.cv)

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
  quant.Yt <- quant.Y[setdiff(1:n.cv,i)]
  quant.Ye <- quant.Y[[i]]
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

  source('_crossvalModel12.R')


  #Compute coverage:
  if(use.cov){
    Sigma.post <- solve(Q.post)
    quant.Yspatial.var <- quant.Xe^2*diag(Sigma.post) + diag(Sigma.E)
  } else {
    quant.Yspatial.var <- quant.Xe^2*diag(inla.qinv(Q.post)) +
              diag(inla.qinv(tau*(obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)))
  }
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
  print(mse)

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


    dev.new()
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
  }
}

