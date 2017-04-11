##
#   comparision method ignoring spatial dependence,
#   fitting the models using qmap
##
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
library(qmap)
#############################
## OPTIONS FOR ESTIMATION
#############################
transfun =  "spline" #"expasympt" #  #"expasympt" #"power.x0", "linear"

model =  
do.plot = TRUE            #Visualise results?
use_BCM_train  = FALSE    #Train model on BCM? If false, use ERA40
use_BCM_eval   = FALSE     #Evaluate results on BCM? If false, use ERA40
n.cv = 8                  #Number of cross validation sets
q = 0.95                  #The quantile to do prediction for
season = 4                #Season to work with (1=winter, 2=spring,...)



source('_data_building_comparisonModel.R')

cat(paste('transform function : ',transfun,'\n',sep=''))

if(transfun == "spline"){ 
  qstep= 0.01 # seems invariant under values of qstep
}else if(transfun == "quant"){
  qstep=0.01 # not as invariant
}else{ 
  qstep=0.02 # seems invariant under values of qstep
}
quant.BCM <- quant.ERA <- quant.Y <- list()
quant.Yp <- list()
model.fit <- list()
##########################
###       TRANING      ##
##########################

for(i in 1:n.cv){
  ind.train <- setdiff(ind,ind.sub[[i]])
  if(use_BCM_train == TRUE){
    train_model <- BCM[,ind.train]
  }else{
    train_model <- ERA[,ind.train]
  }
  Yi <- Y[, ind.train]
  model.fit[[i]] <- list()
  for( j in 1:dim(train_model)[1]) # looping over the grid
  {
   
    if(transfun == "spline"){
      model.fit[[i]][[j]] <- fitQmapSSPLIN(Yi[j, ], train_model[j, ],
                                           wet.day=TRUE,
                                           qstep=qstep)  
    }else if(transfun == "quant"){
      model.fit[[i]][[j]] <- fitQmapQUANT(Yi[j, ], train_model[j, ],
                                           wet.day=TRUE,
                                           qstep=qstep)  
      
    }else{
      model.fit[[i]][[j]] <- fitQmapPTF(Yi[j, ], train_model[j, ],
                                   transfun= transfun,
                                   cost="MAE",wet.day=TRUE,
                                   qstep=qstep)      
    }                    
  }
}

##########################
###       EVALUTION     ##
##########################
model.diff <- list()
MSE <- rep(Inf,n.cv)
for(i in 1:n.cv){
  
  ind.test <- ind.sub[[i]]
  if(use_BCM_eval == TRUE){
    eval_model <- apply(BCM[,ind.test],1,quantile,probs=c(q))
  }else{
    eval_model <- apply(ERA[,ind.test],1,quantile,probs=c(q))
  }
  Yi <- apply(Y[,ind.test],1,quantile,probs=c(q)) 
  model.diff[[i]] <- Inf * rep(1, length(eval_model))
  for( j in 1:length(eval_model)) # looping over the grid
  {
    if(transfun == "spline"){
      Yp <- doQmapSSPLIN(eval_model[j], model.fit[[i]][[j]])
    }else if(transfun == "quant"){
      Yp <- doQmapQUANT(eval_model[j], model.fit[[i]][[j]], ,type="linear")
    }else{
      Yp <- doQmapPTF(eval_model[j], model.fit[[i]][[j]])
    }
    model.diff[[i]][j] <- Yi[j] - Yp                        
  }
  MSE[i] <-  mean(model.diff[[i]]^2)
  cat("    compare model: ", MSE[i], "\n")
}
cat(" total mean : ", mean(MSE), '\n')
# transfun="power.x0", qstep = 0.005
# compare model:  366.116 
# compare model:  293.5718 
# compare model:  235.5294 
# compare model:  207.749 
# compare model:  203.6929 
# compare model:  379.2114 
# compare model:  237.27 
# compare model:  298.6037 
