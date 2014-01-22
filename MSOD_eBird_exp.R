# eBird experiment for predicting the observations (Y)
#
# Author: Jun Yu
# Version: Sep 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSOD.R")
source("MSODL1.R")
source("MSOD_metric.R")

# commandline args
args <- commandArgs(trailingOnly = TRUE)
caseStudy <- args[1]
model <- args[2]
index <- as.numeric(args[3])
if (!(caseStudy %in% c("Hawk","Woodpecker","Finch"))) { stop("Case study species is invalid...\n") }
if (!(model %in% c("OD","ODLP","MSOD"))) { stop("model is invalid...\n") }
if (index < 1 | index > 50) { stop("Index is invalid...\n") }
    
######################
# experiment settings
######################
dataFile <- "./data/eBird_2010_msod.RData"
nRandomRestarts <- 5
if (caseStudy == "Hawk") {
	studySpecies <- c("Accipiter_striatus","Accipiter_cooperii")	
    leakProbs <- c(0.0007,0.01275)
    
    nSpecies <- length(studySpecies)
    confusion <- diag(nSpecies)
    confusion[1,2] <- 1
    confusion[2,1] <- 1
} else if (caseStudy == "Woodpecker") {
	studySpecies <- c("Picoides_pubescens","Picoides_villosus")	
    leakProbs <- c(0.01112,0.003854)
    
    nSpecies <- length(studySpecies)
    confusion <- diag(nSpecies)
    confusion[2,1] <- 1
} else {
	studySpecies <- c("Carpodacus_mexicanus","Carpodacus_purpureus")
    leakProbs <- c(0.0847,0.006141)
    
    nSpecies <- length(studySpecies)
    confusion <- diag(nSpecies)
    confusion[2,1] <- 1
} 

# print the experiment settings
cat("=== Experiment Settings ===\n")
cat("case study:",caseStudy,"\n")
cat("model:",model,"\n")
cat("index:",index,"\n")

#################
# regularization
#################
regType <- 1 # regularization types: 0 for none, 1 for L1, 2 for L2
lambdaO <- 10  # regularization paramters
lambdaD <- 10
lambdaC <-  rep(1,length(studySpecies))
lambda <- list()
lambda$O <- lambdaO
lambda$D <- lambdaD
lambda$C <- lambdaC
cat("=== lambda ===\n")
print(lambda)

#################
# load data
#################
load(dataFile)
dataset <- datasets[[index]]

speciesIdx <- NULL
for (species in studySpecies) {
	speciesIdx <- c(speciesIdx, which(species == speciesNames))
}
print(speciesIdx)

######################
# train and test data
######################
samples <- 1:dataset$nTrSites
trXo     <- dataset$trXo[samples,]
trXd     <- dataset$trXd[samples,,]
trY      <- dataset$trY[samples,,speciesIdx]
trVisits <- dataset$trVisits[samples]
nTrSites <- length(samples) # dataset$nTrSites

teXo     <- dataset$teXo
teXd     <- dataset$teXd
teY      <- dataset$teY[,,speciesIdx]
teVisits <- dataset$teVisits
nTeSites <- dataset$nTeSites

nVisits <- 10

############
# learning
############
teYHat <- array(0, dim = c(nTeSites,nVisits,nSpecies))
if (model == "OD") {
    nRandomRestarts <- 2

#    params <- list()
#    for (s in 1:nSpecies) {
#        
#        # train
#        Y <- trY[,,s]
#        dim(Y) <- c(nrow(Y), ncol(Y), 1)
#        params[[s]] <- RandomRestartEM.MSOD(diag(1),Y,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,c(0)) 
#		
#        # predict        
#        teYHat[,,s] <- PredictDet.MSOD(params[[s]]$params,diag(1),teXo,teXd,teVisits) 
#    } # s
	
	# train
	params <- RandomRestartEM.MSOD(diag(nSpecies),trY,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,rep(0,nSpecies)) 

	# predict
	teYHat <- PredictDet.MSOD(params$params,diag(nSpecies),teXo,teXd,teVisits) 
	
    metric <- evaluateDet(teY,teYHat,teVisits,"True")                    
    print(metric)
} else if (model == "ODLP") {
    nRandomRestarts <- 2

#    params <- list()
#    for (s in 1:nSpecies) {
#        
#        # train
#        Y <- trY[,,s]
#        dim(Y) <- c(nrow(Y), ncol(Y), 1)
#        #params[[s]] <- RandomRestartEM.L1.MSOD(Y,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,NULL) 
#		params[[s]] <- RandomRestartEM.MSOD(diag(1),Y,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,NULL) 
#        
#        # predict        
#        #teYHat[,,s] <- PredictDet.L1.MSOD.R(params[[s]]$params,diag(1),teXo,teXd,teVisits) 
#		teYHat[,,s] <- PredictDet.MSOD(params[[s]]$params,diag(1),teXo,teXd,teVisits) 
#    } # s  
	
	# train
	params <- RandomRestartEM.MSOD(diag(nSpecies),trY,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,NULL) 

	# predict
	teYHat <- PredictDet.MSOD(params$params,diag(nSpecies),teXo,teXd,teVisits) 
	
	
    metric <- evaluateDet(teY,teYHat,teVisits,"True")                    
    print(metric)
} else {    
    # train
    params <- RandomRestartEM.MSOD(confusion,trY,trXo,trXd,trVisits,regType,lambda,nRandomRestarts,NULL) 
    
    # predict
    teYHat <- PredictDet.MSOD(params$params,confusion,teXo,teXd,teVisits) 
}

##############
# save result
##############
resultFile <- paste("./results/eBird/",caseStudy,"_",model,"_",index,".RData",sep="")
save(caseStudy,model,index,teYHat,teY,teVisits,params,file=resultFile)









