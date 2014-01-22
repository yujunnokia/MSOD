# Test MSOD model
#
# Author: Jun Yu
# Version: Feb 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSOD.R")
source("MSOD_genData.R")
source("MSOD_metric.R")

# commandline args
#args <- commandArgs(trailingOnly = TRUE)

######################
# experiment settings
######################
nSpecies <- 3  # number of species
SPECIES_OCCUPANCY_PROB <- runif(3,min=0.2,max=0.8) # c(0.2, 0.4, 0.6)
SPECIES_DETECTION_PROB <- runif(3,min=0.2,max=0.8) # c(0.6, 0.8, 0.5)
LEAK_PROB <- c(0.0,0.0,0.0)
trueConfusion <- diag(nSpecies)
trueConfusion[1,2] <- 1
trueConfusion[2,1] <- 1

# synthetic dataset
nTrSites <- 1000  # number of training sites
nTeSites <- 1000  # number of testing sites
nVisits  <- 3  # number of visits to each site
nOccCovs <- 3  # number of occupancy covariates
nDetCovs <- 4  # number of detection covariates
gaussian <- TRUE # if the covariate is generated from gaussian distribution
nRandomRestarts <- 1  # number of random restarts  

# generate occupancy and detection covariates
trCovs <- GenerateCovariates.MSOD(nTrSites,nVisits,nOccCovs,nDetCovs,gaussian)
trOccCovs <- trCovs$occCovs
trDetCovs <- trCovs$detCovs
teCovs <- GenerateCovariates.MSOD(nTeSites,nVisits,nOccCovs,nDetCovs,gaussian)
teOccCovs <- teCovs$occCovs
teDetCovs <- teCovs$detCovs
occCovs <- rbind(trOccCovs,teOccCovs) 
detCovs <- rbind(matrix(trDetCovs,nrow=nTrSites*nVisits),matrix(teDetCovs,nrow=nTrSites*nVisits))

# set model parameters
trueParams <- list()
for (s in 1:nSpecies) {
    param <- list(alpha=NULL,beta0=NULL,beta=NULL)
    
    # set alpha
    alpha <- rnorm(nOccCovs+1,mean=0,sd=1)
    while( abs(mean((Logistic(occCovs %*% alpha))) - SPECIES_OCCUPANCY_PROB[s]) > 0.001 ) {
        alpha <- rnorm(nOccCovs+1,mean=0,sd=1)
    }
    param$alpha <- alpha
    
    # set leak prob
    param$beta0 <- LEAK_PROB[s]
    
    # set beta
    beta <- rnorm(nDetCovs+1,mean=0,sd=1)
    while( abs(mean((Logistic(detCovs %*% beta))) - SPECIES_DETECTION_PROB[s]) > 0.001 ) {
        beta <- rnorm(nDetCovs+1,mean=0,sd=1)
    }
    cat("Species ",s,"has det prob", mean(Logistic(detCovs %*% beta)), "\n")
    
    param$beta <- list()
    for (r in 1:nSpecies) {
        if (trueConfusion[r,s] == 1) {
            param$beta[[r]] <- rnorm(nDetCovs+1,mean=0,sd=1)
            while( mean((Logistic(detCovs %*% param$beta[[r]]))) > SPECIES_DETECTION_PROB[s]/3 ) {
                param$beta[[r]] <- rnorm(nDetCovs+1,mean=0,sd=1)
            }
            if (s != r) {
                cat("False det prob from",r,"to",s,":", mean(Logistic(detCovs %*% param$beta[[r]])), "\n")
            }
        } else {
            param$beta[[r]] <- rep(0, (nDetCovs+1))
        }
    } # r
    param$beta[[s]] <- beta

    # set the parameter for species s
    trueParams[[s]] <- param
} # s
trueParams$confusion <- trueConfusion

########################
# generate testing data
########################
teVisits <- array(0, nTeSites)
for (i in 1:nTeSites) {
    isMultipleVisits <- runif(1) < 0.5
    if (isMultipleVisits == TRUE) {
        teVisits[i] <- round(runif(1, min=2, max=nVisits))
    } else {
        teVisits[i] <- 1
    }
} # i 
teData <- GenerateData.MSOD(trueParams$confusion,teOccCovs,teDetCovs,teVisits,trueParams[1:nSpecies])
teDetHists <- teData$detHists
teTrueOccs <- teData$trueOccs
cat("Test True Occ avg:\n")
print(colMeans(teTrueOccs))

#########################
# generate training data
#########################
trVisits <- array(0, nTrSites)
for (i in 1:nTrSites) {
    isMultipleVisits <- runif(1) < 0.5
    if (isMultipleVisits == TRUE) {
        trVisits[i] <- round(runif(1, min=2, max=nVisits))
    } else {
        trVisits[i] <- 1
    }
} # i
trData <- GenerateData.MSOD(trueParams$confusion,trOccCovs,trDetCovs,trVisits,trueParams[1:nSpecies])
trDetHists <- trData$detHists
trTrueOccs <- trData$trueOccs
cat("Train True Occ avg:\n")
print(colMeans(trTrueOccs))

# compute overlap sites
cat("Overlap sites:",length(which(rowSums(trTrueOccs)==2)),"\n")

# regularization
regType <- 1
lambda <- list()
lambda$O <- 1
lambda$D <- 1

#############
# learn MSOD
#############
inputConfusion <- matrix(1,nSpecies,nSpecies) # trueConfusion
Y <- trDetHists
Xo <- trOccCovs
Xd <- trDetCovs
visits <- trVisits
leakProbs <- LEAK_PROB

MSODModel <- RandomRestartEM.MSOD(inputConfusion,trDetHists,trOccCovs,trDetCovs,trVisits,regType,lambda,nRandomRestarts,LEAK_PROB) 

###########
# evaluate
###########
for (s in 1:nSpecies) {
    cat("alpha:",s,"\n")
    print(MSODModel$params[[s]]$alpha)
    print(trueParams[[s]]$alpha)

    for (r in 1:nSpecies) {
        if (trueParams$confusion[r,s] == 0) next;
            
        cat("beta:",r,"to",s,"\n")
        print(MSODModel$params[[s]]$beta[[r]])
        print(trueParams[[s]]$beta[[r]])
    }

    cat("beta0:",s,"\n")
    print(MSODModel$params[[s]]$beta0)
    print(trueParams[[s]]$beta0)
}

cat("Learned confusion:\n")
print(MSODModel$params$confusion)
cat("Learned structure:\n")
print(MSODModel$structure)
cat("True confusion:\n")
print(trueConfusion)
cat("Occ and Det rates:\n")
print(MSODModel$occRates)
print(MSODModel$detRates)

# prediction
testPrediction <- Predict.MSOD(MSODModel$params,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
truePrediction <- Predict.MSOD(trueParams,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)

occMetrics <- evaluateOcc(teTrueOccs,testPrediction$probOcc,"MSOD",TRUE)
detMetrics <- evaluateOcc(teTrueOccs,truePrediction$probOcc,"TRUE",TRUE)


# determine the final structure
Y <- teDetHists
Xo <- teOccCovs
Xd <- teDetCovs
visits <- teVisits
finalConfusion <- ComputeStr.MSOD(MSODModel,Y,Xo,Xd,visits)

print(finalConfusion)





