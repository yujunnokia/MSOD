# Implementation of MSOD model
#
# Author: Jun Yu
# Version: Feb 2012
###############################################################################

dyn.load("src/MSODL1.so")

source("MSOD.R")

# maximum number of iterations for each BFGS
MAXIT <- 50
EPSILON <- 1e-10
THRESHOLD <- 1e-5

# regularization on the difference between straight and cross edges
PENALTY <- 1
LAMBDAI <- 10
PENALTYALL <- FALSE

# print trace and debug info
DEBUG <- FALSE
PRINTTRACE <- TRUE

#####################
# Utility functions #
#####################

#
# Compute logistic function
#
Logistic <- function(x) 
{
	y <- 1 / ( 1 + exp(-x))    
	return(y)
}


#
# Implement log(exp(a) + exp(b))
#
LogSumExp <- function(a,b) 
{
	c <- -Inf
	if (b < a) {
		c <- a + log(1 + exp(b - a))
	} else {
		if (a == -Inf && b == -Inf) {
			c <- -Inf
		} else {
			c <- b + log(1 + exp(a - b))
		}
	}    
	return(c)
}

#
# Implement log(exp(a) - exp(b))
#
LogDiffExp <- function(a,b) 
{
	c <- -Inf    
	if (b < a) {
		c <- a + log(1 - exp(b - a))
	} else if (a == b) {
		c <- -Inf
	} else {
		warning("LogDiffExp output -inf.\n")
	}  
	return(c)
}

#
# if b is 0 and c is -inf, then do not add to a.
# else add b * c to a
#
Multiply <- function(a,b)
{
	if ((a == 0 && b == -Inf) || (a == -Inf && b == 0)) {
		return(0)
	}
	return(a*b)
}

#
# find the column index in the probExpectedOccs given the Z
#
FindIdx <- function(Z)
{
	return(sum(unlist(sapply(1:length(Z), function(x,y) 2^(length(y)-x)*y[x], y=Z))) + 1)
}

#
# compute the probability of detection
#
ComputeProbDet <- function(params,Xdit,Z,confusion)
{
	probDet <- 1
    probDet <- probDet * (1-params$beta0)
    for (s in 1:nSpecies) {
        probDet <- probDet * (1-Logistic(Xdit %*% params$beta[[s]]))^(Z[s]*confusion[s])
    } # s
    probDet <- 1 - probDet
	
	if (probDet == 0) { probDet <- EPSILON }
	if (probDet == 1) { probDet <- 1 - EPSILON }
	return(probDet)
}

#
# convert params from array to list
#
ParamsArray2List.L1.MSOD <- function(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs=NULL)
{
	if (class(params) == "list") {
		return(params)
	}
    
	paramsList <- list()
    
	allSpecies <- 1:nSpecies
	confusion <- diag(nSpecies)
	
    # set gamma
	count <- 0
	for (s in 1:nSpecies) {
		confusion[allSpecies[-s],s] <- params[(count+1):(count+nSpecies-1)]
		count <- count + (nSpecies - 1)
	}
    
    # set beta0s
    for (s in 1:nSpecies) {
        paramsList[[s]] <- list()
        
		if (fixedLeak == TRUE) {
			paramsList[[s]]$beta0 <- leakProbs[s];
		} else {
			paramsList[[s]]$beta0 <- params[(count+1):(count+1)]
            count <- count + 1
		}
    }
	
    # set alpha and beta
	for (s in 1:nSpecies) {
        # set alpha
		paramsList[[s]]$alpha <- params[(count+1):(count+nOccCovs)]
		count <- count + nOccCovs
		
        # set beta
        paramsList[[s]]$beta <- list()
        for (ss in 1:nSpecies) {
            if (confusion[ss,s] == 0) { 
                paramsList[[s]]$beta[[ss]] <- rep(0,nDetCovs)
            } else {
                paramsList[[s]]$beta[[ss]] <- params[(count+1):(count+nDetCovs)]
                count <- count + nDetCovs
            }
        } # ss	
	} # s
	
    # set confusion matrix
	paramsList$confusion <- confusion
	
	return(paramsList)
}


##################
# shared methods #
##################

#
# compute marginal log-likelihood for MSOD model
#
ComputeLL.L1.MSOD <- function(params,confusion,Y,Xo,Xd,visits,regType,lambda,idenPenalty=TRUE) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
    lambdaO <- lambda$O
    lambdaD <- lambda$D
    if (is.matrix(lambda$C)) {
        lambdaC <- lambda$C
    } else {
        YY <- NULL
        for (i in 1:nrow(Xo)) {
            tmpY <- Y[i,1:visits[i],1:nSpecies]
            dim(tmpY) <- c(visits[i], nSpecies)
            YY <- rbind(YY, tmpY)
        }
        
        lambdaConfusionMatrix <- matrix(0, nrow=nSpecies, ncol=nSpecies)
        for (s in 1:nSpecies) {
            for (ss in 1:nSpecies) {
                if (s == ss) { next }
                
                lambdaConfusionMatrix[ss,s] <- sum(YY[,ss]==1 & YY[,s]==1) / sum(YY[,ss]==1)
            }
        }
        lambdaC <- lambdaConfusionMatrix * lambda$C
    }
    
	Zs <- expand.grid(rep(list(0:1), nSpecies))

	LL <- 0
	for (i in 1:nSites) {
		Yi    <- matrix(Y[i,,],nrow=nVisits)
		Xoi   <- Xo[i,]
		Xdi   <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		
		iLL <- ComputeiLL.L1.MSOD.C(params,confusion,Yi,Xoi,Xdi,visit,Zs) 	
        #iLL <- ComputeiLL.L1.MSOD.R(params,confusion,Yi,Xoi,Xdi,visit,Zs) 	

		if (iLL == Inf || iLL == -Inf) { stop("iLL is Inf.\n") }
		if (is.na(iLL)) { stop("iLL is NA.\n") }
		if (is.nan(iLL)) { stop("iLL is NaN.\n") }
		
		LL <- LL + iLL
	} # i
    	
	# regularization
	if (regType > 0) {
		for (s in 1:nSpecies) {
            # regularization on alpha
			LL <- LL - 0.5*lambdaO*sum(params[[s]]$alpha[2:nOccCovs]^2) 
			
            ## regularization on beta0
            #LL <- LL - lambdaD*params[[s]]$beta0
            
            # regularization on beta and gamma
            for (ss in 1:nSpecies) {
                # regularization on gamma
                if (ss != s) { 
                    LL <- LL - lambdaC[ss,s]*(confusion[ss,s]^2+EPSILON)^0.5
                }
                
                # regularization on beta
                LL <- LL - 0.5*lambdaD*sum(params[[s]]$beta[[ss]][2:nDetCovs]^2) * confusion[ss,s]
            } # ss
		} # s
	}
	
    if (PENALTY == 1 & idenPenalty) {
        # penalty
        Xdd <- NULL
        for (i in 1:nSites) { Xdd <- rbind(Xdd, Xd[i,1:visits[i],]) }
        
        for (s in 1:nSpecies) {
            avgTrueDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[s]]))
            for (ss in 1:nSpecies) {
                if (ss == s) { next }
                avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[ss]]))
                if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
                    LL <- LL - LAMBDAI * (avgFalseDetProb - avgTrueDetProb)
                }
            }
        }
    }
    
	negLL <- -LL
	if (is.na(LL)) { stop("LL is na...\n") }
    
	return(as.numeric(negLL))
}

#
# compute prob of expected site occupancy for MSOD model
#
ComputeExpectedOccs.L1.MSOD <- function(params,Y,Xo,Xd,visits)
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
    # convert params from array to list
	if (class(params) != "list") {
		stop("parmas must be a list...")
	}
	
	Zs <- expand.grid(rep(list(0:1),nSpecies))
	nZs <- nrow(Zs)
	
	probExpectedOccs <- matrix(0,nrow=nSites,ncol=nZs)
	for (i in 1:nSites) {
		Yi  <- matrix(Y[i,,],nrow=nVisits)
		Xoi <- matrix(Xo[i,],nrow=1)
		Xdi <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		
		denom <- 0
		for (k in 1:nrow(Zs)) {
			Zstr <- paste(Zs[k,],sep="",collapse="")
			Zidx <- FindIdx(Zs[k,])
			
			probExpectedOccs[i,Zidx] <- exp(ComputeiLL.L1.MSOD.C(params[1:nSpecies],params$confusion,Yi,Xoi,Xdi,visit,data.frame(Zs[k,]))) 
			denom <- denom + probExpectedOccs[i,Zidx]
			
			if (is.nan(probExpectedOccs[i,Zidx])) {
				stop("Site",i,"is nan.\n")
			}
		} # k
        if (denom == 0) {
            denom <- 1e-5
        }
		probExpectedOccs[i,] <- probExpectedOccs[i,] / denom
	} # i 
    
    probExpectedOccs[probExpectedOccs == 0] <- 1e-5
    probExpectedOccs[probExpectedOccs == 1] <- 1-1e-5
    
	return(probExpectedOccs)
}


#
# Expectation-Maximization for MSOD model
#
EM.L1.MSOD <- function(Y,Xo,Xd,visits,regType,lambda,initialParams,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
    if (is.null(leakProbs)) {
        fixedLeak <- FALSE
    } else {
        fixedLeak <- TRUE
    }
#	fixedLeak <- FALSE

    if (fixedLeak) {
        params <- c(rep(0.01,nSpecies*(nSpecies-1)),initialParams)   
    } else {
        params <- c(rep(0.01,nSpecies*nSpecies),initialParams)   
    }
    
#	params <- c(rep(1e-2,nSpecies*(nSpecies-1)),initialParams) # assign the initial parameters
#	params <- c(rep(1,nSpecies*(nSpecies-1)),initialParams) # assign the initial parameters 
#	params <- c(runif(nSpecies*(nSpecies-1)),initialParams) # assign the initial parameters
	
	if (PRINTTRACE) { cat("Initial params:",params,"\n") }
	
    # convert params from vector to list
    paramsList <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs)
        
	# E-step
	probExpectedOccs <- ComputeExpectedOccs.L1.MSOD(paramsList,Y,Xo,Xd,visits)
	    
	# EJLL of the initial parameters
	initialEJLL <- -ComputeEJLL.L1.MSOD.C(params,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs)
	newEJLL <- initialEJLL
	newParams <- params
	diffParams <- 1.0e10
	maxIterations <- 50
	iteration <- 1
	tolerance <- 1.0e-5 # 1.0e-5 #0.01
	
	while (diffParams > tolerance && iteration <= maxIterations) {
        # M-step
        if (fixedLeak == FALSE) {
            outputs <- optim(params,ComputeEJLL.L1.MSOD.C,ComputeDerivsOfEJLL.L1.MSOD.C,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs,
            method="L-BFGS-B",
            lower=c(rep(0,nSpecies*nSpecies),rep(-Inf,length(initialParams))),
            upper=c(rep(1,nSpecies*nSpecies),rep(+Inf,length(initialParams))),
            control=list(maxit=MAXIT)) # BFGS  CG  trace=6 in control            
        } else {
            outputs <- optim(params,ComputeEJLL.L1.MSOD.C,ComputeDerivsOfEJLL.L1.MSOD.C,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs,
            method="L-BFGS-B",
            lower=c(rep(0,nSpecies*(nSpecies-1)),rep(-Inf,length(initialParams))),
            upper=c(rep(1,nSpecies*(nSpecies-1)),rep(+Inf,length(initialParams))),
            #				lower=c(rep(0,nSpecies*(nSpecies-1)),rep(-1,length(initialParams))),
            #				upper=c(rep(1,nSpecies*(nSpecies-1)),rep(+1,length(initialParams))),
            #				lower=c(rep(0,nSpecies*(nSpecies-1))),
            #				upper=c(rep(1,nSpecies*(nSpecies-1))),
            control=list(maxit=MAXIT)) # BFGS  CG  trace=6 in control
        }
		params <- outputs$par
		
        # convert params from vector to list
        paramsList <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs)
        
		# E-step
		probExpectedOccs <- ComputeExpectedOccs.L1.MSOD(paramsList,Y,Xo,Xd,visits)
		        
		# udpate params
		oldParams <- newParams
		newParams <- params
		diffParams <- sum((newParams-oldParams)^2) / length(newParams)
		
		# update LL
		oldEJLL <- newEJLL
		newEJLL <- -ComputeEJLL.L1.MSOD.C(newParams,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs)
		
		if (iteration %% 5 == 0 && PRINTTRACE) { 
			newLL <- -ComputeLL.L1.MSOD(paramsList,paramsList$confusion,Y,Xo,Xd,visits,regType,lambda)
			cat("EM.MSOD iteration: ", iteration, " EJLL is ", newEJLL, " LL is ", newLL, "params change is ", diffParams, "\n") 
		}
		iteration <- iteration + 1
	}
	
	return(list(params=paramsList,probExpectedOccs=probExpectedOccs,EJLL=newEJLL))
}

#
# random restart EM for MSOD model
#
RandomRestartEM.L1.MSOD <- function(Y,Xo,Xd,visits,regType=2,lambda,nRandomRestarts=2,leakProbs=NULL,uniformPrior=FALSE) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    YY <- NULL
    Xdd <- NULL
    for (i in 1:nrow(Xo)) {
        tmpY <- Y[i,1:visits[i],1:nSpecies]
        dim(tmpY) <- c(visits[i], nSpecies)
        YY <- rbind(YY, tmpY)
        Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
    }

    lambdaConfusionMatrix <- matrix(0, nrow=nSpecies, ncol=nSpecies)
    for (s in 1:nSpecies) {
        for (ss in 1:nSpecies) {
            if (s == ss) { next }
            cat(sum(YY[,ss]==1 & YY[,s]==1),"/",sum(YY[,ss]==1),"\n")
            if (uniformPrior) {
                lambdaConfusionMatrix[ss,s] <- 1
            } else {
                lambdaConfusionMatrix[ss,s] <- sum(YY[,ss]==1 & YY[,s]==1) / sum(YY[,ss]==1)   
            }
        }
    }
	lambda$C <- lambdaConfusionMatrix * lambda$C
    print(lambda$C)
	    
	# start with full confusion matrix
	confusion <- matrix(1,nrow=nSpecies,ncol=nSpecies)
	
    # if leakPorbs is not specified, then set model to be Noisy-OR-Leak and learn the leak probs
    # o.w. use manually specified leak probs
    if (is.null(leakProbs)) {
        fixedLeak <- FALSE
        cat("Leak probs are NOT fixed...\n")        
    } else {
        fixedLeak <- TRUE
        cat("Leak probs are fixed...\n")
        print(leakProbs)
    }
	
	## fixedLeak is always false and we set leakProbs as the cap in the constraints
	#fixedLeak <- FALSE	
	#if (is.null(leakProbs)) {
#		leakProbs <- rep(1,nSpecies)
	#}
	#cat("Leak probs cap:\n")
	#print(leakProbs)
    
	# determine the number of parameters to be updated 	
	nParams <- nOccCovs*nSpecies + nDetCovs*sum(confusion)
    
	# in iteration 0, we initialize parameters to be the parameters from individual models
	cat("---------------------------\n")
	cat("RandomRestartEM Iteration 0\n")
	initialParams <- array(0.01, c(nParams, 1))
    
	done <- FALSE
	while (!done) 
	{		
		r <- try(bestResult <- EM.L1.MSOD(Y,Xo,Xd,visits,regType,lambda,initialParams,leakProbs))
		done <- !inherits(r, "try-error")
	    initialParams <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1))
	}
	bestEJLL <- bestResult$EJLL
	cat("New EJLL is",bestEJLL,"\n")
	bestLL <- -ComputeLL.L1.MSOD(bestResult$params,bestResult$params$confusion,Y,Xo,Xd,visits,regType,lambda)
	cat("New LL is", bestLL, "\n")
	print(bestResult$params)
	
	# random restart
	for (i in 1:nRandomRestarts) {
		cat("---------------------------\n")
		cat("RandomRestartEML Iteration",i,"\n")
		
		done <- FALSE
		while (!done) 
		{
            initialParams <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1))			
			r <- try(result <- EM.L1.MSOD(Y,Xo,Xd,visits,regType,lambda,initialParams,leakProbs))
			done <- !inherits(r, "try-error")
		}
		
		newEJLL <- result$EJLL 
		cat("Final EJLL is", newEJLL, "\n")
		newLL <- -ComputeLL.L1.MSOD(result$params,result$params$confusion,Y,Xo,Xd,visits,regType,lambda)
		cat("Final LL is", newLL, "\n")
		
        # debug
		print(result$params)
		
#		if (newEJLL > bestEJLL) {
		if (newLL > bestLL) {
			bestResult <- result
			bestLL <- newLL
			bestEJLL <- newEJLL
		}
	}
	
	cat("**********************\n")
	cat("Best EJLL is", bestEJLL, "\n")
    print(bestResult$params)
	
	# compute the marginal log-likelihood
	LL <- ComputeLL.L1.MSOD(bestResult$params,bestResult$params$confusion,Y,Xo,Xd,visits,regType,lambda)
    
    occRates <- NULL
    detRates <- NULL
    for (s in 1:nSpecies) {
        occRates <- c(occRates, mean(Logistic(Xo %*% bestResult$params[[s]]$alpha)))
        detRates <- c(detRates, mean(Logistic(Xdd %*% bestResult$params[[s]]$beta[[s]])))
    }
	
	return(list(params=bestResult$params,confusion=bestResult$confusion,EJLL=bestResult$EJLL,LL=LL,fixedLeak=fixedLeak,
                occRates=occRates, detRates=detRates, obsRates=colMeans(YY)))
}


#
# search for the threshold and recover the structure
#
FindStructure.L1.MSOD <- function(params,Y,Xo,Xd,visits,regType=1,lambda)
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs

	# learned model
	learnedParams <- params$params
	learnedConfusion <- params$params$confusion
	
	YY <- NULL
    Xdd <- NULL
    for (i in 1:nrow(Xo)) {
        tmpY <- Y[i,1:visits[i],1:nSpecies]
        dim(tmpY) <- c(visits[i], nSpecies)
        YY <- rbind(YY, tmpY)
        Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
    }
		
	uniformPrior <- FALSE
	lambdaConfusionMatrix <- matrix(0, nrow=nSpecies, ncol=nSpecies)
	for (s in 1:nSpecies) {
		for (ss in 1:nSpecies) {
			if (s == ss) { next }
			#cat(sum(YY[,ss]==1 & YY[,s]==1),"/",sum(YY[,ss]==1),"\n")
			if (uniformPrior) {
				lambdaConfusionMatrix[ss,s] <- 1
			} else {
				lambdaConfusionMatrix[ss,s] <- sum(YY[,ss]==1 & YY[,s]==1) / sum(YY[,ss]==1)   
			}
		}
	}
	lambda$C <- lambdaConfusionMatrix * lambda$C

	# sort all the pairs in decreasing order
	candidates <- NULL
	for (s in 1:nSpecies) {
		for (ss in 1:nSpecies) {
			if (s == ss) { next }
			candidates <- rbind(candidates, c(s,ss,learnedConfusion[s,ss]))
		}
	}
	candidates <- candidates[order(candidates[,3],decreasing=TRUE),]
	#print(candidates)
	
	# start with empty structure
	confusion <- diag(nSpecies)
	bestLL <- -ComputeLL.L1.MSOD(learnedParams,confusion,Y,Xo,Xd,visits,regType,lambda,FALSE)
	cat("LL:",bestLL,"\n")
	for (i in 1:nrow(candidates)) {
		s <- candidates[i,1]
		ss <- candidates[i,2]
		prob <- candidates[i,3]
        
        if (prob < 0.01) { break }
        
		confusion[s,ss] <- prob
		
		# if the improve is less than some thrshold, then stop
		curLL <- -ComputeLL.L1.MSOD(learnedParams,confusion,Y,Xo,Xd,visits,regType,lambda,FALSE)
		cat(s,"->",ss,"LL:",curLL,"\n")		
		
		if (curLL > bestLL) { 
			bestLL <- curLL
		} else {
			confusion[s,ss] <- 0
			break 
		}
	}
	confusion[confusion > 0] <- 1
	print(confusion)
	
	return(confusion)
}


#
# predict site occupancy
#
PredictOcc.L1.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
{	
	nSpecies <- nrow(confusion)   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	# convert argument params into list
	if (class(params) != "list") {
		stop("parmas must be a list...")
	}
	
	Zs <- expand.grid(rep(list(0:1), nSpecies))
	
	# compute expected occupancy of species s
	probExpectedOccs <- matrix(0,nrow=nSites,ncol=nSpecies)
	for (i in 1:nSites) {
		Yi <- matrix(Y[i,,],nrow=nVisits)
		Xoi <- Xo[i,]
		Xdi <- matrix(Xd[i,,],nrow=nVisits)
		
		denom <- exp(ComputeiLL.L1.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],Zs)) 
		for (s in 1:nSpecies) {
			Z <- Zs[which(Zs[,s]==1),]
			probExpectedOccs[i,s] <- exp(ComputeiLL.L1.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],data.frame(Z))) / denom
					
		} # s
	} # i 
	
	return(probExpectedOccs)
}

PredictDet.L1.MSOD.R <- function(params,confusion,Xo,Xd,visits) 
{
    nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	nSpecies <- nrow(confusion)

    predictedDet <- array(0,c(nSites,nVisits,nSpecies))
    for (s in 1:nSpecies) {
        Zs <- expand.grid(rep(list(0:1), nSpecies))
		nZs <- nrow(Zs)
        
        for (i in 1:nSites) {
            for (t in 1:visits[i]) {
                for (k in 1:nZs) {
                    Z <- Zs[k,]
                    # compute prob of occ
                    probZOcc <- 1
                    for (ss in 1:nSpecies) {
                        probOcc <- Logistic(Xo[i,] %*% params[[ss]]$alpha)
                        probZocc <- probZOcc * (probOcc^Z[ss]) * ((1-probOcc)^(1-Z[ss]))
                    }
                    
                    # comptue prob of det
                    probZDet <- 1 - params[[s]]$beta0
                    for (ss in 1:nSpecies) {
                        probZDet <- probZDet * (1-Logistic(Xd[i,t,] %*% params[[s]]$beta[[ss]]))^(Z[ss]*confusion[ss,s])
                    }
                    probZdet <- 1 - probZDet
                    
                    predictedDet[i,t,s] <- predictedDet[i,t,s] + probZDet*probZOcc
                }
                
                
            }
        }
    }
    
    return(predictedDet)
}

PredictDet.L1.MSOD <- function(params,confusion,Xo,Xd,visits) 
{
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	nSpecies <- nrow(confusion)
	
	# process the parameters
	if (class(params) != "list") {
		stop("parmas must be a list...")
	}
	
	# transpose the matrices
	Xo <- t(Xo)
	Xd <- aperm(Xd,c(3,2,1))
	
	predictedDet <- array(0,c(nSites,nVisits,nSpecies))
	for (s in 1:nSpecies) {
		# extract the params
		alphas <- betas <- beta0s <- NULL
        for (ss in 1:nSpecies) {
            alphas <- c(alphas,params[[ss]]$alpha)
        }
        
        beta0s <- c(beta0s,params[[s]]$beta0)
        
        for (ss in 1:nSpecies) {
            betas <- c(betas,params[[s]]$beta[[ss]])
        }
		
		# all possible parent combinations
		Zs <- expand.grid(rep(list(0:1), nSpecies))
		nZs <- nrow(Zs)
		Zs <- t(Zs)
		
        # make prediction of detection
        result <- .C("PredictDetNoisyORL1", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), 
                     as.integer(nSpecies), as.double(alphas), as.double(betas), as.double(beta0s), as.double(confusion[,s]), as.double(Xo), 
                     as.double(Xd), as.integer(visits), as.integer(Zs), as.integer(nZs),
                     predictedDet = double(nSites*nVisits))
		
		predictedDet[,,s] <- t(array(result[["predictedDet"]],c(nVisits,nSites)))
	} # s
	
	return(predictedDet)
}


#
# predict occ and det
#
Predict.L1.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
{
	probOcc <- PredictOcc.L1.MSOD(params,confusion,Y,Xo,Xd,visits)
	probDet <- PredictDet.L1.MSOD(params,confusion,Xo,Xd,visits)
	
	return(list(probOcc=probOcc,probDet=probDet))
}

#####################
# R and C functions #
#####################

#
# compute marginal log-likelihood for site i for MSOD
# params has to be a list
#
ComputeiLL.L1.MSOD.R <- function(params,confusion,Yi,Xoi,Xdi,visit,Zs) 
{
    nSpecies <- dim(confusion)[1]
	
	if (class(params) != "list") {
		stop("Params has to be a list.")
	}
	
	# compute Occ Probability for each species
	probOccs <- rep(0, nSpecies)
	for (s in 1:nSpecies) {
		probOccs[s] <- Logistic(Xoi %*% params[[s]]$alpha)
	} # s
	
	iL <- 0
	for (k in 1:nrow(Zs)) {
		Z <- Zs[k,]
		isL <- 1
		for (s in 1:nSpecies) {
			isL <- isL * probOccs[s]^Z[s] * (1-probOccs[s])^(1-Z[s])
			
			for ( t in 1:visit) {
				probDet <- ComputeProbDet(params[[s]],Xdi[t,],Z,confusion[,s])
				isL <- isL * probDet^Yi[t,s] * (1-probDet)^(1-Yi[t,s]) 
			} # t
		} # s
		
		iL <- iL + isL
	} # k
	
	iLL <- log(as.numeric(iL))
	return(iLL)
}

ComputeiLL.L1.MSOD.C <- function(params,confusion,Yi,Xoi,Xdi,visit,Zs) 
{
	nOccCovs <- length(Xoi)
	nDetCovs <- dim(Xdi)[2]
	nSpecies <- dim(confusion)[1]
	nZs 	 <- nrow(Zs)
	
	if (class(params) != "list") {
		stop("Params has to be a list.")
	}
	
	# convert params from list to vector
	alphas <- betas <- beta0s <- leakProbs <- NULL
	for (s in 1:nSpecies) {
        alphas <- c(alphas,params[[s]]$alpha)
        
        beta0s <- c(beta0s, params[[s]]$beta0)
        
        for (ss in 1:nSpecies) {
            betas <- c(betas,params[[s]]$beta[[ss]])
        } # ss
    } # s
	
	# transpose Yi and Xdi and Zs
	tYi  <- t(Yi)
	tXdi <- t(Xdi)
	tZs  <- t(Zs)
	
	# call the C subroutine to compute log-likelihood of site i
	result <- .C("ComputeiLLNoisyORL1", as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
                 as.double(alphas), as.double(betas), as.double(beta0s), as.double(confusion), as.integer(tYi),
                 as.double(Xoi), as.double(tXdi), as.integer(visit), as.integer(tZs), as.integer(nZs), 
                 iLL = double(1))
	
	iLL <- result[["iLL"]]
	return(iLL)
}


#
# compute expected joint log-likelihood for MSOD
# the params has to be array
#
ComputeEJLL.L1.MSOD.R <- function(params,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    lambdaO <- lambda$O
    lambdaD <- lambda$D
    lambdaC <- lambda$C
    
	Zs <- expand.grid(rep(list(0:1), nSpecies))
	
	if (class(params) != "list") {
		params <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	confusion <- params$confusion
	
	# compute site occupancy probability for each species
	probOccs <- array(0, c(nSites, nSpecies))
	for (s in 1:nSpecies) {
		probOccs[,s] <- Logistic(Xo %*% params[[s]]$alpha)
	} # s
	
	EJLL <- 0
	for (i in 1:nSites) { 	
		Xdi <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		
		iEJLL <- 0
		for (k in 1:nrow(Zs)) {
			Z <- Zs[k,]
			Zstr <- paste(Zs[k,],sep="",collapse="")
			Zidx <- FindIdx(Z)
			probZ <- probExpectedOccs[i,Zidx]
			
			if (probZ == 0) { next }
			
			isEJLL <- 0
			for (s in 1:nSpecies) {
				isEJLL <- isEJLL + Multiply(Z[s], log(probOccs[i,s])) + Multiply((1-Z[s]), log(1-probOccs[i,s]))
				
				for ( t in 1:visit) {
					probDet <- ComputeProbDet(params[[s]],Xdi[t,],Z,confusion[,s])
					isEJLL <- isEJLL + Multiply(Y[i,t,s], log(probDet)) + Multiply((1-Y[i,t,s]), log(1-probDet)) 
				} # t
			} # s
			iEJLL <- iEJLL + Multiply(probZ,isEJLL)
		} # k
		iEJLL <- as.numeric(iEJLL)
		
		if (iEJLL == Inf || iEJLL == -Inf) { stop("Site ",i," iEJLL is Inf.\n") }
		if (is.na(iEJLL)) { stop("Site ",i," iEJLL is NA.\n") }
		if (is.nan(iEJLL)) { stop("Site ",i," iEJLL is NaN.\n") }
		
		EJLL <- EJLL + iEJLL
	} # i
			
	# regularization
	if (regType > 0) {
		for (s in 1:nSpecies) {
			EJLL <- EJLL - 0.5*lambdaO*sum(params[[s]]$alpha[2:nOccCovs]^2) 
			EJLL <- EJLL - lambdaD*params[[s]]$beta0
            
            for (ss in 1:nSpecies) {
                if (ss != s) {
                    EJLL <- EJLL - lambdaC[ss,s]*(confusion[ss,s]^2 + EPSILON)^0.5
                    #EJLL <- EJLL - 0.5*lambdaC[ss,s]*confusion[ss,s]^2
                }
                
                if (confusion[ss,s] == 0) { next }
                EJLL <- EJLL - 0.5*lambdaD*sum(params[[s]]$beta[[ss]][2:nDetCovs]^2) * confusion[ss,s]
            } # ss            
		} # s
	}
    
	# penalty
    if (PENALTY == 1) {
        Xdd <- NULL
        for (i in 1:nSites) {
            Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
        }
        for (s in 1:nSpecies) {
            avgTrueDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[s]]))
            for (ss in 1:nSpecies) {
                if (ss == s) { next }
                avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[ss]]))
                if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
                    EJLL <- EJLL - LAMBDAI * (avgFalseDetProb - avgTrueDetProb)
                }
            }
        }
    }
    
    #print(EJLL)
    
	negEJLL <- -EJLL
    if (EJLL == Inf || EJLL == -Inf) { stop("EJLL is Inf.\n") }
	if (is.na(EJLL)) { stop("EJLL is NA...\n") }
	if (is.nan(EJLL)) { stop("EJLL is NaN...\n") }    
}

ComputeEJLL.L1.MSOD.C <- function(params,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]  
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
    
    lambdaO <- lambda$O
    lambdaD <- lambda$D
    lambdaC <- lambda$C
	
	if (class(params) != "list") {
		params <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	confusion <- params$confusion
	
	Zs <- expand.grid(rep(list(0:1), nSpecies))
	nZs <- nrow(Zs)
	
	alphas <- betas <- beta0s <- NULL
	for (s in 1:nSpecies) {
        alphas <- c(alphas,params[[s]]$alpha)
        
        beta0s <- c(beta0s, params[[s]]$beta0)
        
        for (ss in 1:nSpecies) {
            betas <- c(betas,params[[s]]$beta[[ss]])
        } # ss
    } # s
	
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tZs <- t(Zs)
	
	# call the C subroutine to compute expected joint log-likelihood
    result <- .C("ComputeEJLLNoisyORL1", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
                 as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.double(confusion), 
                 as.integer(tY), as.double(tXo), as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), 
                 EJLL = double(1))
	EJLL <- result[["EJLL"]]

	# regularization
	if (regType > 0) {
		for (s in 1:nSpecies) {
			EJLL <- EJLL - 0.5*lambdaO*sum(params[[s]]$alpha[2:nOccCovs]^2) 
			EJLL <- EJLL - lambdaD*params[[s]]$beta0
            
            for (ss in 1:nSpecies) {
                if (ss != s) {
                    EJLL <- EJLL - lambdaC[ss,s]*(confusion[ss,s]^2 + EPSILON)^0.5
                    #EJLL <- EJLL - 0.5*lambdaC[ss,s]*confusion[ss,s]^2
                }
                
                if (confusion[ss,s] == 0) { next }
                EJLL <- EJLL - 0.5*lambdaD*sum(params[[s]]$beta[[ss]][2:nDetCovs]^2) * confusion[ss,s]
            } # ss            
		} # s
	}
    	
	# penalty
    if (PENALTY == 1) {
        Xdd <- NULL
        for (i in 1:nSites) {
            Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
        }
        for (s in 1:nSpecies) {
            avgTrueDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[s]]))
            for (ss in 1:nSpecies) {
                if (ss == s) { next }
                avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[ss]]))
                if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
                    EJLL <- EJLL - LAMBDAI * (avgFalseDetProb - avgTrueDetProb)
                }
            }
        }
    }
    
    if (DEBUG) { print(EJLL) }

    if (EJLL == Inf || EJLL == -Inf) { cat("EJLL is Inf.\n") }
	if (is.na(EJLL)) { cat("EJLL is NA...\n") }
	if (is.nan(EJLL)) { cat("EJLL is NaN...\n") }    
    
    negEJLL <- -EJLL
    if (is.infinite(negEJLL) || is.na(negEJLL) || is.nan(negEJLL)) { negEJLL <- 1/EPSILON }
	return(negEJLL)
}


#
# compute derivatives of parameter w.r.t. EJLL for MSOD
# the params has to be an array
#
ComputeDerivsOfEJLL.L1.MSOD.R <- function(params,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    lambdaO <- lambda$O
    lambdaD <- lambda$D
    lambdaC <- lambda$C
    
	Zs <- expand.grid(rep(list(0:1), nSpecies))
    
	if (class(params) != "list") {
		params <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	confusion <- params$confusion
	
	# compute site occupancy probability for each species
	probOccs <- array(0, c(nSites, nSpecies))
	for (s in 1:nSpecies) {
		probOccs[,s] <- Logistic(Xo %*% params[[s]]$alpha)
	}
	
	dQdalpha <- list()
	dQdbeta0 <- list()
	dQdbeta <- list()
	dQdgamma <- matrix(0,nrow=nSpecies,ncol=nSpecies)
	for (s in 1:nSpecies) {
		dQdalpha[[s]] <- rep(0,nOccCovs)
		
        dQdbeta0[[s]] <- 0
        
        dQdbeta[[s]] <- list()
        for (ss in 1:nSpecies) {
            dQdbeta[[s]][[ss]] <- rep(0,nDetCovs)
        } # ss
	} # s
	
	# compute derivatives
	for (s in 1:nSpecies) {
		# derivs of alpha
		for (i in 1:nSites) {
			probExpectedOccis <- 0
			for (idx in which(Zs[,s]==1)) {
				Z <- Zs[idx,]
				Zidx <- FindIdx(Z)
				probExpectedOccis <- probExpectedOccis + probExpectedOccs[i,Zidx]
			} # idx
			
			dQdalpha[[s]] <- dQdalpha[[s]] + (probExpectedOccis-probOccs[i,s])*Xo[i,]
		} # i
		
        # derivs of beta
        for (ss in 1:nSpecies) {
            if (confusion[ss,s] == 0) { next }
            for (i in 1:nSites) {
                for (k in 1:nrow(Zs)) {
                    Z <- Zs[k,]
                    Zidx <- FindIdx(Z)
                    probZ <- probExpectedOccs[i,Zidx]
                    
                    # skip if Z_ss is 0
                    if (Z[ss] == 0) { next }
                    
                    term <- 0
                    for (t in 1:visits[i]) {
                        probDet <- ComputeProbDet(params[[s]],Xd[i,t,],Z,confusion[,s])
                        
                        term <- term + (Y[i,t,s] / probDet - 1) * Logistic(Xd[i,t,] %*% params[[s]]$beta[[ss]]) * Xd[i,t,] * confusion[ss,s]
                    } # t
                    
                    dQdbeta[[s]][[ss]] <- dQdbeta[[s]][[ss]] + probZ * term
                } 		
            } # i
        } # ss
        
        # derivs of gamma
        for (ss in 1:nSpecies) {
            if (ss == s) { next }
            if (confusion[ss,s] < EPSILON) { next }
            
            for (i in 1:nSites) {
                for (k in 1:nrow(Zs)) {
                    Z <- Zs[k,]
                    Zidx <- FindIdx(Z)
                    probZ <- probExpectedOccs[i,Zidx]
                    
                    # skip if Z_ss is 0
                    if (Z[ss] == 0) { next }
                    
                    term <- 0
                    for (t in 1:visits[i]) {
                        probDet <- ComputeProbDet(params[[s]],Xd[i,t,],Z,confusion[,s])
                        term <- term + (1 - Y[i,t,s] / probDet) * log(1-Logistic(Xd[i,t,] %*% params[[s]]$beta[[ss]]))
                    } # t
                    
                    dQdgamma[ss,s] <- dQdgamma[ss,s] + probZ * term
                } # k			
            } # i
        } # ss
        
        # derivs of beta0
		if (fixedLeak == FALSE) {
			for (i in 1:nSites) {
				for (k in 1:nrow(Zs)) {
					Z <- Zs[k,]
					Zidx <- FindIdx(Z)
					probZ <- probExpectedOccs[i,Zidx]
                    
                    term <- 0
                    for (t in 1:visits[i]) {
                        probDet <- ComputeProbDet(params[[s]],Xd[i,t,],Z,confusion[,s])
                        term <- term + (Y[i,t,s] / probDet - 1)
                    } # t
                    
					dQdbeta0[[s]] <- dQdbeta0[[s]] + probZ * term / (1-params[[s]]$beta0)
				} # k			
			} # i
		} # if
        
	} # s
	
	if (regType > 0) {
		for (s in 1:nSpecies) {
			dQdalpha[[s]] <- dQdalpha[[s]] - lambdaO * c(0, params[[s]]$alpha[2:nOccCovs])
			
            if (fixedLeak == FALSE) {
                dQdbeta0[[s]] <- dQdbeta0[[s]] - lambdaD*params[[s]]$beta0
            }
            
            for (ss in 1:nSpecies) {
                dQdbeta[[s]][[ss]]  <- dQdbeta[[s]][[ss]] - lambdaD * c(0, params[[s]]$beta[[ss]][2:nDetCovs])
                dQdgamma[ss,s] <- dQdgamma[ss,s] - lambdaC[ss,s] * (confusion[ss,s]/(sqrt(confusion[ss,s]^2+EPSILON)))
            } # ss
		} # s
	}
	
	# extract the derivs
	derivs <- NULL
	for (s in 1:nSpecies) {
		for (ss in 1:nSpecies) {
			if (ss == s) { next }
			derivs <- c(derivs, dQdgamma[ss,s])
		}
	}

    # beta0
    if (fixedLeak == FALSE) {
        for (s in 1:nSpecies) {
            derivs <- c(derivs, dQdbeta0[[s]])
        }
    }
    
	for (s in 1:nSpecies) {
		# alpha
		derivs <- c(derivs, dQdalpha[[s]])
		
		# beta
        for (ss in 1:nSpecies) {	
            derivs <- c(derivs, dQdbeta[[s]][[ss]])
        } # ss
	} # s
    
	derivs <- -derivs
	return(derivs)
}

ComputeDerivsOfEJLL.L1.MSOD.C <- function(params,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	nSpecies <- dim(Y)[3]
    
    lambdaO <- lambda$O
    lambdaD <- lambda$D
    lambdaC <- lambda$C
	
	if (class(params) != "list") {
		params <- ParamsArray2List.L1.MSOD(params,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	confusion <- params$confusion
	
	Zs <- expand.grid(rep(list(0:1), nSpecies))
	nZs <- nrow(Zs)
	
	alphas <- betas <- beta0s <- NULL
	for (s in 1:nSpecies) {
        alphas <- c(alphas,params[[s]]$alpha)
        
        beta0s <- c(beta0s, params[[s]]$beta0)
        
        for (ss in 1:nSpecies) {
            betas <- c(betas,params[[s]]$beta[[ss]])
        }
    }
	
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tZs <- t(Zs)
	
	# call the C subroutine to compute the derivs of EJLL
	if (fixedLeak == TRUE) {
		result <- .C("ComputeDerivsOfEJLLNoisyORL1", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
				as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.double(confusion), as.integer(tY), as.double(tXo), 
				as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), as.integer(1),
				derivsOfEJLL = double(nSpecies*nOccCovs+nSpecies*nSpecies*nDetCovs+nSpecies*(nSpecies-1)) )
	} else {
		result <- .C("ComputeDerivsOfEJLLNoisyORL1", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
				as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.double(confusion), as.integer(tY), as.double(tXo), 
				as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), as.integer(0),
				derivsOfEJLL = double(nSpecies*(nOccCovs+1)+nSpecies*nSpecies*nDetCovs+nSpecies*(nSpecies-1)) )
	} 
	derivs <- result[["derivsOfEJLL"]]
    
	# penalty
	Xdd <- NULL
	for (i in 1:nSites) {
		Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
	}
	
	# regularization
	if (regType > 0) {
		count <- 1
		for (s in 1:nSpecies) {
			for (ss in 1:nSpecies) {
				if (ss == s) { next }
				derivs[count:count] <- derivs[count:count] - lambdaC[ss,s] * (confusion[ss,s]/sqrt(confusion[ss,s]^2+EPSILON))
				#derivs[count:count] <- derivs[count:count] - lambdaC[ss,s] * confusion[ss,s]
				count <- count + 1
			}
		}
        
        if (fixedLeak==FALSE) {
            for (s in 1:nSpecies) {
                derivs[count:count] <- derivs[count:count] - lambdaD * params[[s]]$beta0
                count <- count + 1
            }
        }
		
		for (s in 1:nSpecies) {
			derivs[(count):(count+nOccCovs-1)] <- derivs[(count):(count+nOccCovs-1)] - lambdaO * c(0, params[[s]]$alpha[2:nOccCovs])
			count <- count + nOccCovs
            
            # find the number of violations on detection probability
            numViolations <- 0
            avgTrueDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[s]]))
            for (ss in 1:nSpecies) {
                avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[ss]]))
                if (avgFalseDetProb > avgTrueDetProb + 1e-10) {
                    numViolations <- numViolations + 1
                }
            }
            if (PENALTYALL) {
                numViolations <- nSpecies-1
            }
            
            for (ss in 1:nSpecies) {
                if (PENALTY == 1) {
                    # penalty
                    avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[ss]]))
                    falseDetProb <- Logistic(Xdd %*% params[[s]]$beta[[ss]])
                    multiplier <- falseDetProb * (1-falseDetProb)
                    if (ss != s) {
                        if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
                            derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] - 
                            LAMBDAI * rowMeans(sapply(1:nrow(Xdd), function(i,a,b) a[i]*b[i,], a = multiplier, b = Xdd))
                        }
                    } else {
                        derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] + 
                        LAMBDAI * rowMeans(sapply(1:nrow(Xdd), function(i,a,b) a[i]*b[i,], a = multiplier, b = Xdd)) * numViolations
                    }
                }
                
                derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] - lambdaD * c(0, params[[s]]$beta[[ss]][2:nDetCovs]) * confusion[ss,s]
                count <- count + nDetCovs
            } # ss
		} # s
	}
	
    if (DEBUG) { print(derivs) }

    if (sum(is.infinite(derivs)) > 0) { print(derivs); stop("derivs contains Inf.\n") }
    if (sum(is.na(derivs)) > 0)       { print(derivs); stop("derivs contains NA...\n") }
    if (sum(is.nan(derivs)) > 0)      { print(derivs); stop("derivs contains NaN...\n") }  

	derivs <- -derivs
	return(derivs)
}
