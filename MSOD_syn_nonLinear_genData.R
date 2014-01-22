# Generate synthetic data from MSOD model but the occupancy status 
# of species can be correlated.
#
# Author: Jun Yu
# Version: Oct 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSODL1.R")
source("MSOD_genData.R")
source("MSOD_metric.R")

######################
# experiment settings
######################
nExps <- 30    # number of experiments

# set number of sites.
nSpecies <- 5    # number of species
nTrSites <- 1000  # number of training sites
nTeSites <- 1000  # number of testing sites
nVisits  <- 3  # number of visits to each site
nOccCovs <- 4  # number of occupancy covariates
nDetCovs <- 4  # number of detection covariates
gaussian <- TRUE # if the covariate is generated from gaussian distribution

datasets <- list()
for (idx in 1:nExps) {
    
    cat("Generate synthetic non-linear dataset",idx,"\n")
    
    ######################
    # generate true model
    ######################
    SPECIES_OCCUPANCY_PROB <- runif(nSpecies,min=0.2, max=0.8) # 0.4
    SPECIES_DETECTION_PROB <- runif(nSpecies,min=0.2, max=0.8) # 0.4
    LEAK_PROB <- rep(0.001,nSpecies)
    trueConfusion <- diag(nSpecies)
    speciesA <- c(sample(1:3,5,replace=TRUE),sample(4:5,3,replace=TRUE))
    speciesB <- c(sample(1:3,5,replace=TRUE),sample(4:5,3,replace=TRUE))
	#speciesA <- sample(1:5,8,replace=TRUE)
	#speciesB <- sample(1:5,8,replace=TRUE)
    for (i in 1:length(speciesA)) {
        trueConfusion[speciesA[i],speciesB[i]] <- 1
    }
    
    #cat("=== True confusion ===\n")
    print(trueConfusion)
    
    ##################################
    # generate occ and det covariates
    ##################################
    trCovs <- GenerateCovariates.MSOD(nTrSites,nVisits,nOccCovs,nDetCovs,gaussian)
    trOccCovs <- trCovs$occCovs
    trDetCovs <- trCovs$detCovs
    teCovs <- GenerateCovariates.MSOD(nTeSites,nVisits,nOccCovs,nDetCovs,gaussian)
    teOccCovs <- teCovs$occCovs
    teDetCovs <- teCovs$detCovs
        
	# make the occupancy nonlinear
	trNLOccCovs <- trOccCovs
	trNLDetCovs <- trDetCovs
	teNLOccCovs <- teOccCovs
	teNLDetCovs <- teDetCovs
	
	trNLOccCovs[,1] <- sin(trOccCovs[,1]+1)
	trNLOccCovs[,2] <- exp(trOccCovs[,2]-1)
	trNLOccCovs[,3] <- trOccCovs[,3]*trOccCovs[,4]

	teNLOccCovs[,1] <- sin(teOccCovs[,1]+1)
	teNLOccCovs[,2] <- exp(teOccCovs[,2]-1)
	teNLOccCovs[,3] <- teOccCovs[,3]*teOccCovs[,4]

    occCovs <- rbind(trNLOccCovs,teNLOccCovs) 
    detCovs <- rbind(matrix(trNLDetCovs,nrow=nTrSites*nVisits),matrix(teNLDetCovs,nrow=nTrSites*nVisits))
	
    # set model parameters
    trueParams <- list()
    for (s in 1:nSpecies) {
        param <- list(alpha=rnorm(nOccCovs+1,mean=0,sd=1),beta0=NULL,beta=NULL)
        
        # set alpha
        while( abs(mean((Logistic(occCovs %*% param$alpha))) - SPECIES_OCCUPANCY_PROB[s]) > 0.001 ) {
            param$alpha <- rnorm(nOccCovs+1,mean=0,sd=1)
        }
        
        # set leak prob
        param$beta0 <- LEAK_PROB[s]
        
        # set beta
        beta <- rnorm(nDetCovs+1,mean=0,sd=1)
        while( abs(mean((Logistic(detCovs %*% beta))) - SPECIES_DETECTION_PROB[s]) > 0.001 ) {
            beta <- rnorm(nDetCovs+1,mean=0,sd=1)
        }
        #cat("Species ",s,"has det prob", mean(Logistic(detCovs %*% beta)), "\n")
        
        param$beta <- list()
        for (ss in 1:nSpecies) {
            if (trueConfusion[ss,s] == 1) {
                param$beta[[ss]] <- rnorm(nDetCovs+1,mean=0,sd=1)
                while( mean((Logistic(detCovs %*% param$beta[[ss]]))) > SPECIES_DETECTION_PROB[s]*0.7 |
                       mean((Logistic(detCovs %*% param$beta[[ss]]))) < SPECIES_DETECTION_PROB[s]*0.3) {
                    param$beta[[ss]] <- rnorm(nDetCovs+1,mean=0,sd=1)
                }
                if (s != ss) {
                    #cat("False det prob from",ss,"to",s,":", mean(Logistic(detCovs %*% param$beta[[ss]])), "\n")
                }
            } else {
                param$beta[[ss]] <- rep(0, (nDetCovs+1))
            }
        } # ss
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
    teData <- GenerateData.MSOD(trueParams$confusion,teNLOccCovs,teNLDetCovs,teVisits,trueParams[1:nSpecies])
    teDetHists <- teData$detHists
    teTrueOccs <- teData$trueOccs
    #cat("Test True Occ avg:\n")
    #print(colMeans(teTrueOccs))
    
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
    trData <- GenerateData.MSOD(trueParams$confusion,trNLOccCovs,trNLDetCovs,trVisits,trueParams[1:nSpecies])
    trDetHists <- trData$detHists
    trTrueOccs <- trData$trueOccs
    #cat("Train True Occ avg:\n")
    #print(colMeans(trTrueOccs))
    
    dataset <- list(index=idx,nSpecies=nSpecies,trueParams=trueParams,gaussian=gaussian,nVisits=nVisits,
                    SPECIES_OCCUPANCY_PROB=SPECIES_OCCUPANCY_PROB, SPECIES_DETECTION_PROB=SPECIES_DETECTION_PROB, LEAK_PROB=LEAK_PROB,
					trNLOccCovs=trNLOccCovs,teNLOccCovs=teNLOccCovs,trNLDetCovs=trNLDetCovs,teNLDetCovs=teNLDetCovs,
                    nTrSites=nTrSites,trOccCovs=trOccCovs,trDetCovs=trDetCovs,trVisits=trVisits,trDetHists=trDetHists,trTrueOccs=trTrueOccs,
                    nTeSites=nTeSites,teOccCovs=teOccCovs,teDetCovs=teDetCovs,teVisits=teVisits,teDetHists=teDetHists,teTrueOccs=teTrueOccs)
    datasets[[idx]] <- dataset
} # i

# save datasets
save(datasets,file="./data/MSOD_syndata_nonLinear.RData")



