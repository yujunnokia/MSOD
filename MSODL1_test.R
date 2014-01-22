# Test different methods in MSOD model on synthetic data
#
# Author: Jun Yu
# Version: Feb 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/eBird.MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/eBird.MSOD")
#source("MSOD.R")
source("MSODL1.R")
source("MSOD.genData.R")
source("MSOD.metric.R")

## commandline args
#args <- commandArgs(trailingOnly = TRUE)
#occ1 <- as.numeric(args[1])
#occ2 <- as.numeric(args[2])
#trueLeakProb <- as.numeric(args[3])

######################
# experiment settings
######################

# set number of parameters and confusion matrix
nSpecies <- 1
SPECIES_OCCUPANCY_PROB <- c(0.2)
LEAK_PROB <- c(0.2)
trueConfusion <- diag(nSpecies)

#nSpecies <- 2  # number of species
#SPECIES_OCCUPANCY_PROB <- c(0.5, 0.2) 
#LEAK_PROB <- c(0.1,0)
#trueConfusion <- diag(nSpecies)
#trueConfusion[1,2] <- 1
##trueConfusion[2,1] <- 1

#nSpecies <- 3  # number of species
#SPECIES_OCCUPANCY_PROB <- c(0.6,0.3,0.1)
#LEAK_PROB <- c(0.2,0.01,0.01)
#trueConfusion <- diag(nSpecies)
##trueConfusion[1,2] <- 1
#trueConfusion[2,3] <- 1
##trueConfusion[3,2] <- 1

#nSpecies <- 4  # number of species
#SPECIES_OCCUPANCY_PROB <- c(0.6, 0.5, 0.3, 0.2) 

#nSpecies <- 5  # number of species
#SPECIES_OCCUPANCY_PROB <- c(0.6, 0.5, 0.4, 0.3, 0.2) 
#trueConfusion <- diag(nSpecies)
#trueConfusion[1,2] <- 1
#trueConfusion[2,3] <- 1
##trueConfusion[3,1] <- 1
##trueConfusion[3,2] <- 1
#trueConfusion[4,5] <- 1
##trueConfusion[1,5] <- 1
##trueConfusion[5,3] <- 1

#nSpecies <- 3  # number of species
#trueConfusion <- diag(nSpecies)
#for (i in 1:floor(1 * nSpecies)) {
#	s1 <- sample(1:nSpecies, 1)
#	s2 <- sample(1:nSpecies, 1)
#	trueConfusion[s1,s2] <- 1
#}
cat("=== True confusion ===\n")
print(trueConfusion)

# set number of sites..
nTrSites <- 500  # number of training sites
nTeSites <- 500  # number of testing sites
nVisits  <- 3  # number of visits to each site
nOccCovs <- 3  # number of occupancy covariates
nDetCovs <- 4  # number of detection covariates
gaussian <- TRUE # if the covariate is generated from gaussian distribution
nRandomRestarts <- 2  # number of random restarts  

# generate occupancy and detection covariates
trCovs <- GenerateCovariates.MSOD(nTrSites,nVisits,nOccCovs,nDetCovs,gaussian)
trOccCovs <- trCovs$occCovs
trDetCovs <- trCovs$detCovs
teCovs <- GenerateCovariates.MSOD(nTeSites,nVisits,nOccCovs,nDetCovs,gaussian)
teOccCovs <- teCovs$occCovs
teDetCovs <- teCovs$detCovs
occCovs <- rbind(trOccCovs,teOccCovs) 
detCovs <- rbind(matrix(trDetCovs,nrow=nTrSites*nVisits),matrix(teDetCovs,nrow=nTrSites*nVisits))

# set leakProbs
cat("=== leakProbs ===\n")
print(LEAK_PROB)

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
    while( abs(mean((Logistic(detCovs %*% beta))) - 0.5) > 0.001 ) {
        beta <- rnorm(nDetCovs+1,mean=0,sd=1)
	}
    cat("Species ",s,"has det prob", mean(Logistic(detCovs %*% beta)), "\n")
	
    param$beta <- list()
    for (ss in 1:nSpecies) {
        if (trueConfusion[ss,s] == 1) {
            param$beta[[ss]] <- rnorm(nDetCovs+1,mean=0,sd=1)
            while( mean((Logistic(detCovs %*% param$beta[[ss]]))) > mean((Logistic(detCovs %*% beta))) ) {
                param$beta[[ss]] <- rnorm(nDetCovs+1,mean=0,sd=1)
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

################
# Testing cases
################
Y <- trDetHists
Xd <- trDetCovs
Xo <- trOccCovs
visits <- trVisits

###################
# Test on Syn data
###################

# test whether the model can learn the leak probability from data
if (1 == 1) {    
    # regularization
    regType <- 1 # regularization types: 0 for none, 1 for L1, 2 for L2
    lambdaO <- 1  # regularization paramters
    lambdaD <- 1#10
    lambdaC <-  rep(1,nSpecies) #rep(0.01,nSpecies)
    lambda <- list()
    lambda$O <- lambdaO
    lambda$D <- lambdaD
    lambda$C <- lambdaC
    cat("=== lambda ===\n")
    cat("Alpha:",lambda$O,"Beta:",lambda$D,"Gamma:",lambda$C,"\n")
    
    # learn MSOD model
    learntParams <- RandomRestartEM.L1.MSOD(Y,Xo,Xd,visits,regType,lambda,nRandomRestarts) 
    
    for (s in 1:nSpecies) {
        cat("alpha:",s,"\n")
        print(learntParams$params[[s]]$alpha)
        print(trueParams[[s]]$alpha)
        
        for (ss in 1:nSpecies) {
            if (trueParams$confusion[ss,s] == 0) next;
            cat("beta:",ss,"to",s,"\n")
            print(learntParams$params[[s]]$beta[[ss]])
            print(trueParams[[s]]$beta[[ss]])
        }
        
        cat("beta0:",s,"\n")
        print(learntParams$params[[s]]$beta0)
        print(trueParams[[s]]$beta0)
    }
    
    confusion <- learntParams$params$confusion
    confusion[confusion > 1e-1] <- 1
    confusion[confusion <= 1e-1] <- 0
    cat("Learned confusion:\n")
    print(confusion)
    cat("True confusion:\n")
    print(trueConfusion)
    
    cat("True LL: ", -ComputeLL.L1.MSOD(trueParams,trueConfusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC),"\n") 
    cat("MSOD.L1 LL: ", -ComputeLL.L1.MSOD(learntParams$params,learntParams$params$confusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC),"\n") 
    
    testPrediction <- Predict.L1.MSOD(learntParams$params,learntParams$params$confusion,teDetHists,teOccCovs,teDetCovs,teVisits)
    truePrediction <- Predict.L1.MSOD(trueParams,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
    
    evaluateOcc(teTrueOccs,testPrediction$probOcc,"MSOD")
    evaluateOcc(teTrueOccs,truePrediction$probOcc,"TRUE")
    
    print(lambda)
    cat("Lambda on adding constraints: ",LAMBDAI," \n")
    cat("==================================================\n")
}

# Test EM.L1.MSOD function
if (1 == 1) {    
    # regularization
    regType <- 1 # regularization types: 0 for none, 1 for L1, 2 for L2
    lambdaO <- 1  # regularization paramters
    lambdaD <- 10
    lambdaC <-  rep(1,nSpecies) #rep(0.01,nSpecies)
    lambda <- list()
    lambda$O <- lambdaO
    lambda$D <- lambdaD
    lambda$C <- lambdaC
    cat("=== lambda ===\n")
    cat("Alpha:",lambda$O,"Beta:",lambda$D,"Gamma:",lambda$C,"\n")
    
    # learn MSOD model with different fixed leak probs
    fixedLeakProbs <- c(0.1) # c(0.0,0.1,0.2,0.4)
    for (fixedLeakProb in fixedLeakProbs)
    {    
        testLeakProbs <- rep(0,nSpecies)  # leakProbs
        testLeakProbs[1] <- fixedLeakProb
        
        learntParams <- RandomRestartEM.L1.MSOD(Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,testLeakProbs) 
        
        for (s in 1:nSpecies) {
            cat("alpha:",s,"\n")
            print(learntParams$params[[s]]$alpha)
            print(trueParams[[s]]$alpha)
            
            for (ss in 1:nSpecies) {
                if (trueParams$confusion[ss,s] == 0) next;
                cat("beta:",ss,"to",s,"\n")
                print(learntParams$params[[s]]$beta[[ss]])
                print(trueParams[[s]]$beta[[ss]])
            }
            
            cat("beta0:",s,"\n")
            print(learntParams$params[[s]]$beta0)
            print(trueParams[[s]]$beta0)
        }
        
        confusion <- learntParams$params$confusion
        confusion[confusion > 1e-1] <- 1
        confusion[confusion <= 1e-1] <- 0
        cat("Learned confusion:\n")
        print(confusion)
        cat("True confusion:\n")
        print(trueConfusion)
        
        cat("True LL: ", -ComputeLL.L1.MSOD(trueParams,trueConfusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC),"\n") 
        cat("MSOD.L1 LL: ", -ComputeLL.L1.MSOD(learntParams$params,learntParams$params$confusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC),"\n") 
        
        testPrediction <- Predict.L1.MSOD(learntParams$params,learntParams$params$confusion,teDetHists,teOccCovs,teDetCovs,teVisits)
        truePrediction <- Predict.L1.MSOD(trueParams,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
        
        evaluateOcc(teTrueOccs,testPrediction$probOcc,"MSOD")
        evaluateOcc(teTrueOccs,truePrediction$probOcc,"TRUE")
        
        print(lambda)
        cat("Lambda on adding constraints: ",LAMBDAI," \n")
        cat("==================================================\n")
    }
}


##########################
## Test single function ##
##########################
# Test ComputeiLL.MSOD function
if (1 == 0) {
	i <- 2
	Yi  <- matrix(Y[i,,],nrow=nVisits)
	Xoi <- matrix(Xo[i,],nrow=1)
	Xdi <- matrix(Xd[i,,],nrow=nVisits)
	visit <- visits[i]
	Zs <- matrix(sample(0:1,nSpecies,replace=TRUE),ncol=nSpecies)
	
	confusion <- trueParams$confusion
	confusion[1,2] <- 0.5
	confusion[2,1] <- 0.7
    	
	# params,confusion,Yi,Xoi,Xdi,visit,Zs,model
	cat("ComputeiLL.L1.MSOD.R: ", ComputeiLL.L1.MSOD.R(trueParams[1:nSpecies],confusion,Yi,Xoi,Xdi,visit,Zs),"\n") 
	cat("ComputeiLL.L1.MSOD.C: ", ComputeiLL.L1.MSOD.C(trueParams[1:nSpecies],confusion,Yi,Xoi,Xdi,visit,Zs),"\n") 
}

# Test ComputeLL.MSOD function
if (1 == 0) {
	cat("ComputeLL.L1.MSOD: ", -ComputeLL.L1.MSOD(trueParams[1:nSpecies],trueParams$confusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC),"\n")  
}

# Test ComputeEJLL.MSOD function
if (1 == 0) {
	probExpectedOccs <- ComputeExpectedOccs.L1.MSOD(trueParams,Y,Xo,Xd,visits)
	
    #	trueParams$confusion[1,2] <- 0.82
    #	trueParams$confusion[2,1] <- 0.28
	
	cat("ComputeEJLL.L1.MSOD.R: ", -ComputeEJLL.L1.MSOD.R(trueParams,probExpectedOccs,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC,FALSE),"\n")  
	cat("ComputeEJLL.L1.MSOD.C: ", -ComputeEJLL.L1.MSOD.C(trueParams,probExpectedOccs,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,lambdaC,FALSE),"\n")  
}

# Test ComputeDerivsOfEJLL.MSOD function
if (1 == 0) {
	probExpectedOccs <- ComputeExpectedOccs.L1.MSOD(trueParams,Y,Xo,Xd,visits)
	
	trueParams$confusion[1,2] <- 0.82
	trueParams$confusion[2,1] <- 0.28
	
	cat("ComputeDerivsOfEJLL.L1.MSOD.R: ", -ComputeDerivsOfEJLL.L1.MSOD.R(trueParams,probExpectedOccs,Y,Xo,Xd,visits,0,lambdaO,lambdaD,lambdaC,FALSE),"\n")
	cat("ComputeDerivsOfEJLL.L1.MSOD.C: ", -ComputeDerivsOfEJLL.L1.MSOD.C(trueParams,probExpectedOccs,Y,Xo,Xd,visits,0,lambdaO,lambdaD,lambdaC,FALSE),"\n")
}

## Test RandomRestartEM.MSOD function
#{
#	result <- RandomRestartEM.MSOD(trueConfusion,Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,model,initialParams=NULL)
#	params <- result$params
#	ComputeParamsDiff(trueParams,params,trueConfusion,model) 
#	
#	# output log likelihood
#	cat("==== Log Likelihood ====\n")
#	cat("MSOD LL:",-ComputeLL.MSOD(params,trueConfusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,model),"\n")
#	cat("True LL:",-ComputeLL.MSOD(trueParams,trueConfusion,Y,Xo,Xd,visits,regType,lambdaO,lambdaD,model),"\n")
#	
#	# predict occ and det
#	tePrediction   <- Predict.MSOD(params,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits,model)
#	truePrediction <- Predict.MSOD(trueParams,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits,model)
#	
#	cat("==== Occupancy Evaluation ====\n")
#	evaluateOcc(teTrueOccs,tePrediction$probOcc,"MSOD")
#	evaluateOcc(teTrueOccs,truePrediction$probOcc,"TRUE")
#}
