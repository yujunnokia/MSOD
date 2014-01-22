# Implementation of MSOD model given a fixed confusion matrix
#
# Author: Jun Yu
# Version: June 2013
###############################################################################

dyn.load("src/MSOD.so")

# maximum number of iterations for each BFGS
BFGS_MAX_ITERS <- 100
EPSILON <- 1e-10
MAX_ITERATIONS <- 50
TOLERANCE <- 1.0e-5

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
ParamsArray2List.MSOD <- function(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs=NULL)
{
	if (class(params) == "list") {
		return(params)
	}
    
	paramsList <- list()
	count <- 0
    
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
	
	return(paramsList)
}


##################
# shared methods #
##################

#
# compute marginal log-likelihood for MSOD model
#
ComputeLL.MSOD <- function(params,confusion,Y,Xo,Xd,visits,regType,lambda) 
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

	LL <- 0
	for (i in 1:nSites) {
		Yi    <- matrix(Y[i,,],nrow=nVisits)
		Xoi   <- Xo[i,]
		Xdi   <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		
		iLL <- ComputeiLL.MSOD.C(params,confusion,Yi,Xoi,Xdi,visit,Zs) 	
        #iLL <- ComputeiLL.MSOD.R(params,confusion,Yi,Xoi,Xdi,visit,Zs) 	

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
			
            # regularization on beta0
            LL <- LL - lambdaD*params[[s]]$beta0
            
            # regularization on beta
            for (ss in 1:nSpecies) {
                LL <- LL - 0.5*lambdaD*sum(params[[s]]$beta[[ss]][2:nDetCovs]^2)  * confusion[ss,s]
            } # ss
		} # s
	}
	
    if (PENALTY == 1) {
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
ComputeExpectedOccs.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
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
			
			probExpectedOccs[i,Zidx] <- exp(ComputeiLL.MSOD.C(params[1:nSpecies],confusion,Yi,Xoi,Xdi,visit,data.frame(Zs[k,]))) 
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
    
    probExpectedOccs[probExpectedOccs == 0] <- EPSILON
    probExpectedOccs[probExpectedOccs == 1] <- 1-EPSILON
    
	return(probExpectedOccs)
}


#
# Expectation-Maximization for MSOD model
#
EM.MSOD <- function(confusion,Y,Xo,Xd,visits,regType,lambda,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	nParams <- nOccCovs*nSpecies + nDetCovs*sum(confusion)
    if (is.null(leakProbs)) {
        fixedLeak <- FALSE
        nParams <- nParams + nSpecies
    } else {
        fixedLeak <- TRUE
    }
    
    
    comb <- sapply(0:(2^nSpecies-1),function(x){ as.integer(intToBits(x))})[nSpecies:1,]
    if (nSpecies == 1) {
        dim(comb) <- c(1, 2)
    }

    # initialize Zs
	probExpectedOccs <- array(1,c(nSites,2^nSpecies))
    for (i in 1:nSites) {
    	if (nSpecies == 1) {
            Yi <- c(sum(Y[i,,1]))
    	} else {
            Yi <- colSums(Y[i,,])
    	}
        
        for (s in 1:nSpecies) {
            Z1 <- which(comb[s,] == 1)
            Z0 <- which(comb[s,] == 0)
            
            if (Yi[s] > 0) {
                probExpectedOccs[i,Z1] <- probExpectedOccs[i,Z1] * 0.9
                probExpectedOccs[i,Z0] <- probExpectedOccs[i,Z0] * 0.1               
            } else {
                probExpectedOccs[i,Z1] <- probExpectedOccs[i,Z1] * 0.3
                probExpectedOccs[i,Z0] <- probExpectedOccs[i,Z0] * 0.7
            }
        }
    }
    #print(Y[1:5,,])
    #print(probExpectedOccs[1:5,])
	    
	# EJLL of the initial parameters
	newEJLL <- -Inf
	params <- array(rnorm(nParams,mean=0,sd=1), c(nParams, 1)) #  array(0, c(nParams, 1)) 
    cat("Initial params:",params,"\n")
	diffParams <- 1/EPSILON
	iteration <- 1	
	while (diffParams > TOLERANCE && iteration <= MAX_ITERATIONS) {
        oldParams <- params
        
        # M-step
        if (fixedLeak == FALSE) {
            outputs <- optim(params,ComputeEJLL.MSOD.C,ComputeDerivsOfEJLL.MSOD.C,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs,
            method="L-BFGS-B",
            lower=c(rep(0,nSpecies),rep(-Inf,nParams)),
			#upper=c(leakProbs,rep(+Inf,nParams)),
            upper=c(rep(1,nSpecies),rep(+Inf,nParams)),
            control=list(maxit=BFGS_MAX_ITERS)) # BFGS  CG  trace=6 in control            
        } else {
            outputs <- optim(params,ComputeEJLL.MSOD.C,ComputeDerivsOfEJLL.MSOD.C,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs,
            method="BFGS",
            control=list(maxit=BFGS_MAX_ITERS)) # BFGS  CG  trace=6 in control
        }
		params <- outputs$par
		
        # convert params from vector to list
        paramsList <- ParamsArray2List.MSOD(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs)
        
		# E-step
		probExpectedOccs <- ComputeExpectedOccs.MSOD(paramsList,confusion,Y,Xo,Xd,visits)
		        
		# udpate params
		diffParams <- sum((params-oldParams)^2) / length(params)
				
		if (iteration %% 5 == 0 && PRINTTRACE) { 
			newLL <- -ComputeLL.MSOD(paramsList,confusion,Y,Xo,Xd,visits,regType,lambda)
            newEJLL <- -ComputeEJLL.MSOD.C(paramsList,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs)
			cat("EM.MSOD iteration: ", iteration, " EJLL is ", newEJLL, " LL is ", newLL, "params change is ", diffParams, "\n") 
		}
		iteration <- iteration + 1
	}
	
	return(list(params=paramsList,probExpectedOccs=probExpectedOccs,LL=newLL,EJLL=newEJLL))
}

#
# random restart EM for MSOD model
#
RandomRestartEM.MSOD <- function(confusion,Y,Xo,Xd,visits,regType=1,lambda,nRandomRestarts=2,leakProbs=NULL,uniformPrior=FALSE) 
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
    
	# random restart
    bestLL <- -Inf
	for (i in 1:nRandomRestarts) {
		cat("---------------------------\n")
		cat("RandomRestartEM Iteration",i,"\n")
		
		done <- FALSE
		while (!done) 
		{
			r <- try(result <- EM.MSOD(confusion,Y,Xo,Xd,visits,regType,lambda,leakProbs))
			done <- !inherits(r, "try-error")
		}
		
		newEJLL <- result$EJLL 
		cat("Final EJLL is", newEJLL, "\n")
		newLL <- -ComputeLL.MSOD(result$params,confusion,Y,Xo,Xd,visits,regType,lambda)
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
	    
    occRates <- array(0,c(1,nSpecies))
    detRates <- array(0,c(nSpecies,nSpecies))
    for (s in 1:nSpecies) {
        occRates[s] <- mean(Logistic(Xo %*% bestResult$params[[s]]$alpha))
        for (r in 1:nSpecies) {
            if (confusion[r,s] == 0) {
                detRates[r,s] <- 0
            } else {
                detRates[r,s] <- mean(Logistic(Xdd %*% bestResult$params[[s]]$beta[[r]]))
            }
        }
    }
	
	return(list(params=bestResult$params,confusion=confusion,
                EJLL=bestResult$EJLL,LL=bestResult$LL,
                regType=regType,lambda=lambda,fixedLeak=fixedLeak,
                occRates=occRates, detRates=detRates, obsRates=colMeans(YY)))
}

#
# use a validation set to determine the structure of the MSOD model
#
ComputeStr.MSOD <- function(model,Y,Xo,Xd,visits) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    scores <- data.frame(from=rep(0,nSpecies*(nSpecies-1)),
                         to=rep(0,nSpecies*(nSpecies-1)),
                         detProb=rep(0,nSpecies*(nSpecies-1)))
    n <- 1
    for (r in 1:nSpecies) {
        for (s in 1:nSpecies) {
            if (r == s) { next }
            scores$from[n] <- r
            scores$to[n] <- s
            scores$detProb[n] <- model$detRates[r,s]
            n <- n+1
        }
    }
    
    # sort all the unspecified cross edges by their detection probabilities
    scores <- scores[order(-scores$detProb),,drop=FALSE]
    print(scores)
    
    # and add them in a greedy fashion
    confusion <- diag(nSpecies)
    bestLL <- -ComputeLL.MSOD(model$params,confusion,Y,Xo,Xd,visits,model$regType,model$lambda)
    n <- 1
    while (n <= nrow(scores)) {
        confusion[scores[n,1],scores[n,2]] <- 1
        LL <- -ComputeLL.MSOD(model$params,confusion,Y,Xo,Xd,visits,model$regType,model$lambda)
        cat(LL,"\n")
        if (LL >= bestLL) {
            bestLL <- LL
        } else {
            confusion[scores[n,1],scores[n,2]] <- 0
            break
        }
        n <- n + 1
    }
    
    return(confusion)
}

#
# predict site occupancy
#
PredictOcc.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
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
		
		denom <- exp(ComputeiLL.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],Zs)) 
		for (s in 1:nSpecies) {
			Z <- Zs[which(Zs[,s]==1),]
			probExpectedOccs[i,s] <- exp(ComputeiLL.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],data.frame(Z))) / denom
					
		} # s
	} # i 
	
	return(probExpectedOccs)
}

PredictDet.MSOD <- function(params,confusion,Xo,Xd,visits) 
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
        result <- .C("PredictDetNoisyOR", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), 
                     as.integer(nSpecies), as.double(alphas), as.double(betas), as.double(beta0s), as.integer(confusion[,s]), as.double(Xo), 
                     as.double(Xd), as.integer(visits), as.integer(Zs), as.integer(nZs),
                     predictedDet = double(nSites*nVisits))
		
		predictedDet[,,s] <- t(array(result[["predictedDet"]],c(nVisits,nSites)))
	} # s
	
	return(predictedDet)
}


#
# predict occ and det
#
Predict.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
{
	probOcc <- PredictOcc.MSOD(params,confusion,Y,Xo,Xd,visits)
	probDet <- PredictDet.MSOD(params,confusion,Xo,Xd,visits)
	
	return(list(probOcc=probOcc,probDet=probDet))
}

#####################
# R and C functions #
#####################

#
# compute marginal log-likelihood for site i for MSOD
# params has to be a list
#
ComputeiLL.MSOD.R <- function(params,confusion,Yi,Xoi,Xdi,visit,Zs) 
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

ComputeiLL.MSOD.C <- function(params,confusion,Yi,Xoi,Xdi,visit,Zs) 
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
	result <- .C("ComputeiLLNoisyOR", as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
                 as.double(alphas), as.double(betas), as.double(beta0s), as.integer(confusion), as.integer(tYi),
                 as.double(Xoi), as.double(tXdi), as.integer(visit), as.integer(tZs), as.integer(nZs), 
                 iLL = double(1))
	
	iLL <- result[["iLL"]]
	return(iLL)
}


#
# compute expected joint log-likelihood for MSOD
# the params has to be array
#
ComputeEJLL.MSOD.R <- function(params,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
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
		params <- ParamsArray2List.MSOD(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
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

ComputeEJLL.MSOD.C <- function(params,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
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
		params <- ParamsArray2List.MSOD(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}

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
    result <- .C("ComputeEJLLNoisyOR", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
                 as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.integer(confusion), 
                 as.integer(tY), as.double(tXo), as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), 
                 EJLL = double(1))
	EJLL <- result[["EJLL"]]

	# regularization
	if (regType > 0) {
		for (s in 1:nSpecies) {
			EJLL <- EJLL - 0.5*lambdaO*sum(params[[s]]$alpha[2:nOccCovs]^2) 
			EJLL <- EJLL - lambdaD*params[[s]]$beta0
            
            for (ss in 1:nSpecies) {
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

    #if (EJLL == Inf || EJLL == -Inf) { cat("EJLL is Inf.\n") }
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
ComputeDerivsOfEJLL.MSOD.R <- function(params,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
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
		params <- ParamsArray2List.MSOD(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
	# compute site occupancy probability for each species
	probOccs <- array(0, c(nSites, nSpecies))
	for (s in 1:nSpecies) {
		probOccs[,s] <- Logistic(Xo %*% params[[s]]$alpha)
	}
	
	dQdalpha <- list()
	dQdbeta0 <- list()
	dQdbeta <- list()
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
            } # ss
		} # s
	}
	
	# extract the derivs
	derivs <- NULL

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
            if (confusion[ss,s] == 0) { next }
            derivs <- c(derivs, dQdbeta[[s]][[ss]])
        } # ss
	} # s
    
	derivs <- -derivs
	return(derivs)
}

ComputeDerivsOfEJLL.MSOD.C <- function(params,confusion,probExpectedOccs,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
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
		params <- ParamsArray2List.MSOD(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
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
		result <- .C("ComputeDerivsOfEJLLNoisyOR", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
				as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.integer(confusion), as.integer(tY), as.double(tXo), 
				as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), as.integer(1),
				derivsOfEJLL = double(nSpecies*nOccCovs+sum(confusion)*nDetCovs) )
	} else {
		result <- .C("ComputeDerivsOfEJLLNoisyOR", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
				as.double(alphas), as.double(betas), as.double(beta0s), as.double(probExpectedOccs), as.integer(confusion), as.integer(tY), as.double(tXo), 
				as.double(tXd), as.integer(visits), as.integer(tZs), as.integer(nZs), as.integer(0),
				derivsOfEJLL = double(nSpecies*(nOccCovs+1)+sum(confusion)*nDetCovs) )
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
                if (confusion[ss,s] == 0) { next }

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


########################
# non-independent MSOD
########################

ComputeiLL.NI.MSOD.C <- function(params,confusion,Yi,Xoi,Xdi,visit,Zs) 
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
	result <- .C("ComputeiLLNoisyORNonIndependent", as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
                 as.double(alphas), as.double(betas), as.double(beta0s), as.integer(confusion), as.integer(tYi),
                 as.double(Xoi), as.double(tXdi), as.integer(visit), as.integer(tZs), as.integer(nZs), 
                 iLL = double(1))
	
	iLL <- result[["iLL"]]
	return(iLL)
}


#
# predict site occupancy
#
PredictOcc.NI.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
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
		
		denom <- exp(ComputeiLL.NI.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],Zs)) 
		for (s in 1:nSpecies) {
			Z <- Zs[which(Zs[,s]==1),]
			probExpectedOccs[i,s] <- exp(ComputeiLL.NI.MSOD.C(params,confusion,Yi,Xoi,Xdi,visits[i],data.frame(Z))) / denom
					
		} # s
	} # i 
	
	return(probExpectedOccs)
}

PredictDet.NI.MSOD <- function(params,confusion,Xo,Xd,visits) 
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
        result <- .C("PredictDetNoisyORNonIndependent", as.integer(nSites), as.integer(nVisits), as.integer(nOccCovs), as.integer(nDetCovs), 
                     as.integer(nSpecies), as.double(alphas), as.double(betas), as.double(beta0s), as.integer(confusion[,s]), as.double(Xo), 
                     as.double(Xd), as.integer(visits), as.integer(Zs), as.integer(nZs),
                     predictedDet = double(nSites*nVisits))
		
		predictedDet[,,s] <- t(array(result[["predictedDet"]],c(nVisits,nSites)))
	} # s
	
	return(predictedDet)
}


#
# predict occ and det
#
Predict.NI.MSOD <- function(params,confusion,Y,Xo,Xd,visits)
{
	probOcc <- PredictOcc.NI.MSOD(params,confusion,Y,Xo,Xd,visits)
	probDet <- PredictDet.NI.MSOD(params,confusion,Xo,Xd,visits)
	
	return(list(probOcc=probOcc,probDet=probDet))
}
