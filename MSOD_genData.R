# Generate synthetic data from the Multi-Species OD model
#
# Author: Jun Yu
# Version: Sep 2012
###############################################################################

#source("MSOD.R")


#
# generate occ and det covariates
#
GenerateCovariates.MSOD <- function(nSites,nVisits,nOccCovs,nDetCovs,gaussian=TRUE) 
{
	# generate occupancy and detection covariates
	if (gaussian) {
		occCovs <- rnorm(nSites * (nOccCovs+1))
		detCovs <- rnorm(nSites * nVisits * (nDetCovs+1))
	} else {
		occCovs <- runif(nSites * (nOccCovs+1))*2-1
		detCovs <- runif(nSites * nVisits * (nDetCovs+1))*2-1
	}
	dim(occCovs) <- c(nSites, (nOccCovs+1))
	dim(detCovs) <- c(nSites, nVisits, (nDetCovs+1))
	occCovs[,1] <- array(1, c(nSites, 1))
	detCovs[,,1] <- array(1, c(nSites, nVisits, 1))
	
	return(list(occCovs=occCovs,detCovs=detCovs))
}

#
# generate data of the MSOD model
# 
GenerateData.MSOD <- function(confusion,occCovs,detCovs,visits,params) 
{
	nSpecies <- length(params)
	nSites   <- dim(occCovs)[1]
	nOccCovs <- dim(occCovs)[2]
	nDetCovs <- dim(detCovs)[3]
	nVisits  <- dim(detCovs)[2]
	
	# generate true occupancy
	trueOccs <- matrix(0,nrow=nSites,ncol=nSpecies)
	for (s in 1:nSpecies) {
		trueOccs[,s] <- runif(nSites) < Logistic(occCovs %*% params[[s]]$alpha)
	} # s
	trueOccs[trueOccs == TRUE] <- 1
	
	# generate observations
	detHists <- array(0, c(nSites,nVisits,nSpecies))
	for (i in 1:nSites) {
		Z <- trueOccs[i,]
		
		for (t in 1:visits[i]) {
			for (s in 1:nSpecies) {
				pi <- which(confusion[,s] == 1)
				
				# compute detection prob
				probDet <- ComputeProbDet(params[[s]],detCovs[i,t,],Z,confusion[,s])
				
				# generate observations for species s
				if (runif(1) < probDet) {
					detHists[i,t,s] <- 1
				}	
			} # s		
		} # t
	} # i
	
	return(list(trueOccs=trueOccs, detHists=detHists))
}

#
# generate synthetic data from non-independent MSOD model
# 
GenerateData.nonIndependent.MSOD <- function(confusion,occCovs,detCovs,visits,params,occupancyMatrix) 
{
	nSpecies <- length(params)
	nSites   <- dim(occCovs)[1]
	nOccCovs <- dim(occCovs)[2]
	nDetCovs <- dim(detCovs)[3]
	nVisits  <- dim(detCovs)[2]
	
	# generate true occupancy
	trueOccs <- matrix(0,nrow=nSites,ncol=nSpecies)
	for (s in 1:nSpecies) {
		trueOccs[,s] <- runif(nSites) < Logistic(occCovs %*% params[[s]]$alpha)
	} # s
	trueOccs[trueOccs == TRUE] <- 1
	
	childSpecies <- which(colSums(occupancyMatrix) != 0)
	#print(childSpecies)
	#print(colMeans(trueOccs))
	for (s in childSpecies) {
		competitor <- which(occupancyMatrix[,s] == -1)
		mutualism <- which(occupancyMatrix[,s] == 1)
		
		if (length(competitor) > 0) { 
			cat("competitor",competitor," -> ",s,"\n")
			for (i in 1:nSites) {
				if (trueOccs[i,competitor] == 1) {
					trueOccs[i,s] <- runif(1) < (Logistic(occCovs[i,] %*% params[[s]]$alpha)/2)
				}
			}
		} 
		
		if (length(mutualism) > 0) { 
			cat("mutualism",mutualism," -> ",s,"\n")
			for (i in 1:nSites) {
				if (trueOccs[i,mutualism] == 1) {
					trueOccs[i,s] <- runif(1) < (Logistic(occCovs[i,] %*% params[[s]]$alpha)*1.2)
				}
			}
		} 		
	}
	#print(colMeans(trueOccs))
		
	# generate observations
	detHists <- array(0, c(nSites,nVisits,nSpecies))
	for (i in 1:nSites) {
		Z <- trueOccs[i,]
		
		for (t in 1:visits[i]) {
			for (s in 1:nSpecies) {
				pi <- which(confusion[,s] == 1)
				
				# compute detection prob
				probDet <- ComputeProbDet(params[[s]],detCovs[i,t,],Z,confusion[,s])
				
				# generate observations for species s
				if (runif(1) < probDet) {
					detHists[i,t,s] <- 1
				}	
			} # s		
		} # t
	} # i
	
	return(list(trueOccs=trueOccs, detHists=detHists))
}

