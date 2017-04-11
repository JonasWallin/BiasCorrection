cat("-- Model details --\n")
cat("Model : ")
cat("\nTraining Data :")
if(use_BCM_train){cat(" BCM")}else{cat(" ERA40")}
cat("\nEvaluation Data :")
if(use_BCM_eval){cat(" BCM ")}else {cat(" ERA40 ")}
cat("\n")

#############################
## LOAD DATA
#############################
files <- c("index.norway.1x1kmAggregert.sav",
           "lat.norway.sav",
           "lon.norway.sav",
           "dataTrainingNotRandom.obs.sav",
           "dataTrainingNotRandom.bcm.sav",
           "dataTrainingNotRandom.era.sav",
           "dataTestNotRandom.bcm.sav",
           "dataTestNotRandom.obs.sav",
           "dataTestNotRandom.era.sav")
if(is.null(data_location)){
  data.url <- "http://www.math.chalmers.se/~bodavid/data/downscaling/"
  cat('Downloading files may take sometime\n')
  for(i in 1:length(files))
  {
    cat('*')
    con <- url(paste(data.url,files[i],sep=""))
    load(con); close(con)
  }
  cat('Donwloading files done\n')
}else{
  for(i in 1:length(files))
    load(paste(data_location,files[i],sep=""))
}
m = 58
n = 63
obs.ind = (index.norway.1x1kmAggregert[,2]-1)*m+index.norway.1x1kmAggregert[,1]


n.grid <- length(dataTraining.bcm[[season]])
n.obs.train <- length(dataTraining.bcm[[season]][[1]])
n.obs.eval <- length(dataTest.bcm[[season]][[1]])

Y_train <- ERA_train <- BCM_train <- matrix(0,n.grid,n.obs.train)
Y_eval <- ERA_eval <- BCM_eval <- matrix(0,n.grid,n.obs.eval)


for(t in 1:n.grid)
{
  Y_train[t,] = dataTraining.obs[[season]][[t]]
  Y_eval[t,]  =dataTest.obs[[season]][[t]]

  BCM_train[t,] = dataTraining.bcm[[season]][[t]]
  BCM_eval[t,]  =dataTest.bcm[[season]][[t]]

  ERA_train[t,] = dataTraining.era[[season]][[t]]
  ERA_eval[t,]  =dataTest.era[[season]][[t]]
}

BCM <- cBind(BCM_train,BCM_eval)
Y   <- cBind(Y_train, Y_eval)
ERA <- cBind(ERA_train,ERA_eval)


###################################
## Extract quantiles for each fold
###################################

n <- dim(Y)[2]
ind <- 1:n
n.sub = ceiling(n/n.cv)
ind.sub <- list()

tmp <- matrix(1:(n.sub*n.cv),n.sub,n.cv)
for(i in 1:(n.cv-1)){
  ind.sub[[i]] = tmp[,i]
}
ind.sub[[n.cv]] = intersect(tmp[,n.cv],ind)



