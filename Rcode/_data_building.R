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

#############################
## transform data if selected
#############################
if(use_log)
{
  BCM  = log(BCM)
  Y  = log(Y)
  ERA = log(ERA)
}

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

quant.BCM <- quant.ERA <- quant.Y <- list()
for(i in 1:n.cv){
  quant.BCM[[i]] <- apply(BCM[,ind.sub[[i]]],1,quantile,probs=c(q))
  quant.ERA[[i]] <- apply(ERA[,ind.sub[[i]]],1,quantile,probs=c(q))
  quant.Y[[i]]   <- apply(Y[,ind.sub[[i]]],1,quantile,probs=c(q))
}


#########################################
## Construct matrices for spatial models
#########################################

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

load_smooth <- function(type, data_location = NULL)
{
  #data 
  #type is 'BCM' or 'ERA'
  smooth.X = TRUE
  
  if(type=='BCM'){
    use_BCMs = TRUE
  }else{
    use_BCMs = FALSE
  }

  if(is.null(data_location)==FALSE){
    file_name = paste(data_location,type)
    if(use_log)
      file_name <- paste(file_name,"_log",sep="")
    
    file_name <- paste(file_name,"_season_", season, sep = "")
    file_name <- paste(file_name,".RData", sep="")
    if(file.exists(file_name)){
      load(file_name)
    }else{
      source('smoothing_X.R', local=TRUE)
    }
    
  }else{
    source('smoothing_X.R', local=TRUE)
  }
  #
  return(X_smooth)
}