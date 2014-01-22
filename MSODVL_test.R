# Test MSOD model with variational learning
#
# Author: Jun Yu
# Version: Nov 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSOD.R")
source("MSODVL.R")
source("MSOD_genData.R")
source("MSOD_metric.R")

# commandline args
#args <- commandArgs(trailingOnly = TRUE)

######################
# experiment settings
######################
nSpecies <- 2  # number of species
SPECIES_OCCUPANCY_PROB <- runif(nSpecies,min=0.2,max=0.8) # c(0.2, 0.4, 0.6)
SPECIES_DETECTION_PROB <- runif(nSpecies,min=0.2,max=0.8) # c(0.6, 0.8, 0.5)
LEAK_PROB <- rep(0.0001,nSpecies)
trueConfusion <- diag(nSpecies)
trueConfusion[1,2] <- 1
trueConfusion[2,1] <- 1
#trueConfusion[3,2] <- 1
#trueConfusion[4,5] <- 1

print(trueConfusion)

# synthetic dataset
nTrSites <- 1000  # number of training sites
nTeSites <- 1000  # number of testing sites
nVisits  <- 3    # number of visits to each site
nOccCovs <- 3    # number of occupancy covariates
nDetCovs <- 3  # number of detection covariates
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
            while( mean((Logistic(detCovs %*% param$beta[[r]]))) > SPECIES_DETECTION_PROB[s]/2 ) {
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

################
# Testing cases
################
inputConfusion <- matrix(1,nSpecies,nSpecies) # trueConfusion
Y <- trDetHists
Xd <- trDetCovs
Xo <- trOccCovs
visits <- trVisits
leakProbs <- LEAK_PROB
Z <- trTrueOccs

# Test UpdateVariationalParams 
{
	cat("=== Test UpdateVariationalParams function ===\n")
	params <- trueParams[1:nSpecies]
	confusion <- trueParams$confusion
	Z <- trTrueOccs
	
	q.R <- UpdateVariationalParams.R(params,Z,confusion,Y,Xo,Xd,visits)
    print(q.R[1:5,1,,1])
	q.C <- UpdateVariationalParams.C(params,Z,confusion,Y,Xo,Xd,visits)
    print(q.C[1:5,1,,1])
}

# Test ComputeExpectedOccs.MSODVL
{
	cat("=== Test ComputeExpectedOccs function ===\n")
	params <- trueParams[1:nSpecies]
	confusion <- trueParams$confusion
	Z <- trTrueOccs
	
	q <- UpdateVariationalParams.R(params,Z,confusion,Y,Xo,Xd,visits)
	
	Zhat.R <- ComputeExpectedOccs.MSODVL.R(params,q,confusion,Y,Xo,Xd,visits)
	print(Zhat.R[1:5,])
	Zhat.C <- ComputeExpectedOccs.MSODVL.C(params,q,confusion,Y,Xo,Xd,visits)
	print(Zhat.C[1:5,])
	
	print(Z[1:5,])
}

# Test ComputeEJLL.MSODVL
{
	cat("=== Test ComputeEJLL.MSODVL function ===\n")
	cat("EJLL VL.R:",ComputeEJLL.MSODVL.R(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,leakProbs),"\n")
	cat("EJLL VL.C:",ComputeEJLL.MSODVL.C(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,leakProbs),"\n")
	
	nParams <- (nOccCovs+1)*nSpecies + (nDetCovs+1)*nSpecies*nSpecies
	newParams <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1))
	cat("EJLL VL.R:",ComputeEJLL.MSODVL.R(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,leakProbs),"\n")
	cat("EJLL VL.C:",ComputeEJLL.MSODVL.C(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,leakProbs),"\n")
	
#	probExpectedOccs <- ComputeExpectedOccs.MSOD(params,confusion,Y,Xo,Xd,visits) 
#	cat("EJLL:",ComputeEJLL.MSOD.R(params,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB),"\n")
}

# Test ComputeDerivsOfEJLL.MSODVL
{
	cat("=== Test ComputeDerivsOfEJLL.MSODVL function ===\n")
	derivs.R <- ComputeDerivsOfEJLL.MSODVL.R(params,q,confusion,Z,Y,Xo,Xd,visits,0,lambda,TRUE,leakProbs)
    derivs.C <- ComputeDerivsOfEJLL.MSODVL.C(params,q,confusion,Z,Y,Xo,Xd,visits,0,lambda,TRUE,leakProbs)
    print(derivs.R - derivs.C)

	nParams <- (nOccCovs+1)*nSpecies + (nDetCovs+1)*nSpecies*nSpecies
	newParams <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1))
	derivs.R <- ComputeDerivsOfEJLL.MSODVL.R(newParams,q,confusion,Z,Y,Xo,Xd,visits,0,lambda,TRUE,leakProbs)
	derivs.C <- ComputeDerivsOfEJLL.MSODVL.C(newParams,q,confusion,Z,Y,Xo,Xd,visits,0,lambda,TRUE,leakProbs)
    print(derivs.R - derivs.C)
}

# Test M-step
{
	cat("=== Test M-step in EM ===\n")
	# E-step: first update q and then calcualte expected Z
	q <- UpdateVariationalParams.R(params,Z,confusion,Y,Xo,Xd,visits)
	Zhat <- ComputeExpectedOccs.MSODVL.C(params,q,confusion,Y,Xo,Xd,visits)
	
##	cat("EJLL VL.R:",ComputeEJLL.MSODVL.R(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB),"\n")
##	print(ComputeDerivsOfEJLL.MSODVL.R(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB))
##	
##	cat("EJLL VL.C:",ComputeEJLL.MSODVL.C(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB),"\n")
##	print(ComputeDerivsOfEJLL.MSODVL.C(newParams,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB))
	
	
	# M-step: update the model parameters
#	outputs <- optim(newParams,ComputeEJLL.MSODVL.R,ComputeDerivsOfEJLL.MSODVL.R,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB,
#			method="BFGS",
#			control=list(maxit=MAXIT)) # BFGS  CG  trace=6 in control
#	print(ParamsArray2List.MSODVL(outputs$par,confusion,nSpecies,nOccCovs+1,nDetCovs+1,TRUE,LEAK_PROB))
	
    nParams <- (nOccCovs+1)*nSpecies + (nDetCovs+1)*nSpecies*nSpecies
    newParams <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1))
	outputs <- optim(newParams,ComputeEJLL.MSODVL.C,ComputeDerivsOfEJLL.MSODVL.C,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,TRUE,LEAK_PROB,
			method="BFGS",
			control=list(maxit=BFGS_MAX_ITERS)) # BFGS  CG  trace=6 in control
    learnedParams <- ParamsArray2List.MSODVL(outputs$par,confusion,nSpecies,nOccCovs+1,nDetCovs+1,TRUE,LEAK_PROB)
	print(learnedParams[[1]]$alpha)
	print(params[[1]]$alpha)
	print(learnedParams[[2]]$alpha)
	print(params[[2]]$alpha)
}


{
    MSOD <- RandomRestartEM.MSOD(inputConfusion,Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,leakProbs)
    MSODVL.true <- RandomRestartEM.MSODVL(inputConfusion,Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,leakProbs,trTrueOccs)
    MSODVL <- RandomRestartEM.MSODVL(inputConfusion,Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,leakProbs)
    
    # prediction
    MSODPrediction <- Predict.MSOD(MSOD$params,inputConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
    MSODVLPrediction.true <- Predict.MSODVL(MSODVL.true$params,inputConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
    MSODVLPrediction <- Predict.MSODVL(MSODVL$params,inputConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)
    truePrediction <- Predict.MSOD(trueParams,trueConfusion,teDetHists,teOccCovs,teDetCovs,teVisits)

    evaluateOcc(teTrueOccs,MSODPrediction$probOcc,"MSOD",TRUE)
    evaluateOcc(teTrueOccs,MSODVLPrediction.true$probOcc,"MSODVL.true",TRUE)
    evaluateOcc(teTrueOccs,MSODVLPrediction$probOcc,"MSODVL",TRUE)
    evaluateOcc(teTrueOccs,truePrediction$probOcc,"TRUE",TRUE)
    
	evaluateDet(teDetHists,MSODPrediction$probDet,teVisits,"MSOD",TRUE)
    evaluateDet(teDetHists,MSODVLPrediction.true$probDet,teVisits,"MSODVL.true",TRUE)
    evaluateDet(teDetHists,MSODVLPrediction$probDet,teVisits,"MSODVL",TRUE)
    evaluateDet(teDetHists,truePrediction$probDet,teVisits,"TRUE",TRUE)
}

if (TRUE == FALSE) {
	cat("=== Test EM ===\n")
	confusion <- matrix(1,nSpecies,nSpecies)
	for (restart in 1:1) {
		result <- EM.MSODVL(confusion,Y,Xo,Xd,visits,regType,lambda,LEAK_PROB)
		cat("--- Not give Z ---\n")
		print(result$params[[1]]$alpha)
		print(trueParams[[1]]$alpha)
		print(result$params[[2]]$alpha)
		print(trueParams[[2]]$alpha)
	}
	
	evaluateOcc(Z,result$Z,"MSODVL",TRUE)
	
	truePrediction  <- PredictOcc.MSOD(trueParams,trueConfusion,Y,Xo,Xd,visits)
	evaluateOcc(Z,truePrediction,"true",TRUE)
	
	result <- EM.MSODVL(confusion,Y,Xo,Xd,visits,regType,lambda,LEAK_PROB,Z)
	cat("--- Give Z ---\n")
	print(result$params[[1]]$alpha)
	print(trueParams[[1]]$alpha)
	print(result$params[[2]]$alpha)
	print(trueParams[[2]]$alpha)
	
	Y <- teDetHists
	Xd <- teDetCovs
	Xo <- teOccCovs
	visits <- teVisits
	
	# predict Z
	modelPrediction <- PredictOcc.MSODVL(result$params,confusion,Y,Xo,Xd,visits)
	truePrediction  <- PredictOcc.MSOD(trueParams,trueConfusion,Y,Xo,Xd,visits)
	
	# evaluate Z
	occMetrics <- evaluateOcc(teTrueOccs,modelPrediction,"MSODVL",TRUE)
	occMetrics <- evaluateOcc(teTrueOccs,truePrediction,"TRUE",TRUE)
}


#{
#	x = 1
#	k = 0
#	while (k < 100) {
#		x = 1 + 1/x
#		print(x)
#		k = k+1
#	}
#}
