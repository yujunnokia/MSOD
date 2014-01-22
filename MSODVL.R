# Implementation of MSODVL model given a fixed confusion matrix
#
# Author: Jun Yu
# Version: June 2013
###############################################################################

dyn.load("src/MSODVL.so")

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
# f function 
#
F <- function(x)
{
	return(log(1-exp(-x)))
}


#
# compute the probability of detection
#
ComputeProbDet.MSODVL <- function(params,q,Xdit,Z,confusion)
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
ParamsArray2List.MSODVL <- function(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs=NULL)
{
	if (class(params) == "list") {
		return(params)
	}
	
	paramsList <- list()
	
	# set beta0s
	for (r in 1:nSpecies) {
		paramsList[[r]] <- list()
		paramsList[[r]]$beta <- list()
		paramsList[[r]]$beta0 <- leakProbs[r];
	} # r
	
	# set alpha and beta
	count <- 0
	for (s in 1:nSpecies) {
		# set alpha
		paramsList[[s]]$alpha <- params[(count+1):(count+nOccCovs)]
		count <- count + nOccCovs
		
		# set beta
		for (r in 1:nSpecies) {
			paramsList[[s]]$beta[[r]] <- params[(count+1):(count+nDetCovs)]
			count <- count + nDetCovs
		} # ss    
	} # s
	
	return(paramsList)
}


##################
# shared methods #
##################

#
# compute marginal log-likelihood for MSODVL model
#
ComputeLL.MSODVL <- function(params,q,confusion,Y,Xo,Xd,visits,regType,lambda) 
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
		
		iLL <- ComputeiLL.MSODVL.C(params,confusion,Yi,Xoi,Xdi,visit,Zs)     
		#iLL <- ComputeiLL.MSODVL.R(params,confusion,Yi,Xoi,Xdi,visit,Zs)     
		
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
# compute prob of expected site occupancy for MSODVL model
# return a list of Z1 and Z0
# Z1: P(Z_ir=1 | X_i) \prod_t \prod_s P(Y_its | Z_ir=1, W_it)
# Z0: P(Z_ir=0 | X_i) \prod_t \prod_s P(Y_its | Z_ir=0, W_it)
#
ComputeExpectedOccs.MSODVL.R <- function(params,q,confusion,Y,Xo,Xd,visits)
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	# convert params from array to list
	if (class(params) != "list") { stop("parmas must be a list...") }
	
	Z <- matrix(0,nrow=nSites,ncol=nSpecies)
	for (i in 1:nSites) {
		Yi  <- matrix(Y[i,,],nrow=nVisits)
		Xoi <- matrix(Xo[i,],nrow=1)
		Xdi <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		qi <- q[i,,,]
		
		Z[i,] <- ComputeiLL.MSODVL.R(params,qi,confusion,Yi,Xoi,Xdi,visit)
	} # i 
	
	Z[Z == 0] <- EPSILON
	Z[Z == 1] <- 1-EPSILON
	
	return(Z)
}

ComputeExpectedOccs.MSODVL.C <- function(params,q,confusion,Y,Xo,Xd,visits)
{
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	nSpecies <- dim(Y)[3]
	
	if (class(params) != "list") { stop("params must be a list...") }
	
	alphas <- betas <- beta0s <- NULL
	for (r in 1:nSpecies) {
		alphas <- c(alphas,params[[r]]$alpha)
		
		beta0s <- c(beta0s, params[[r]]$beta0)
		
		for (s in 1:nSpecies) {
			betas <- c(betas,params[[r]]$beta[[s]])
		}
	}
	
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tq  <- aperm(q,c(4,3,2,1))
	
	# call the C subroutine to compute the derivs of EJLL
	result <- .C("ComputeExpectedOccMSODVL", 
			as.integer(nSites), 
			as.integer(nVisits), 
			as.integer(nOccCovs), 
			as.integer(nDetCovs), 
			as.integer(nSpecies),
			as.double(alphas), 
			as.double(betas), 
			as.double(beta0s), 
			as.double(tq), 
			as.integer(tY), 
			as.double(tXo), 
			as.double(tXd), 
			as.integer(visits),
			Z = double(nSites*nSpecies) )
	Z <- result[["Z"]]
	
#	print(Z)
	
	# reshape
	dim(Z) <- c(nSpecies,nSites)
	Z <- aperm(Z,c(2,1))
	
	return(Z)
}


#
# update variational parameter q
# use fixed-point iterative method
#
UpdateVariationalParams.R <- function(params,Z,confusion,Y,Xo,Xd,visits) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	if (class(params) != "list") { stop("parmas must be a list...") }
	
	theta0 <- rep(0,nSpecies)
	for (s in 1:nSpecies) {
		theta0[s] <- -log(1-params[[s]]$beta0)
        if (theta0[s] == 0) { theta0[s] <- EPSILON }
	}
	
	q <- array(0,c(nSites,nVisits,nSpecies,nSpecies))
	for (i in 1:nSites) {
		for (t in 1:visits[i]) {
			for (s in 1:nSpecies) {
				for (r in 1:nSpecies) {
                    ditrs <- Logistic(Xd[i,t,] %*% params[[s]]$beta[[r]])
                    if (ditrs == 1) { ditrs <- 1-EPSILON }
                    if (ditrs == 0) { ditrs <- EPSILON }
					theta  <- -log(1-ditrs)
                    
                    #cat("i:",i,"t:",t,"s:",s,"r:",r,"d:",ditrs,"\n")
					
					qitrs <- 0.5
					k <- 0
					while (k < 50) {
						A <- exp(-theta0[s]-theta/qitrs)
                        if (A == 0) { A <- EPSILON}
                        if (A == 1) { A <- 1-EPSILON}
                        F0 <- F(theta0[s])
                            						
						#if (i == 1) { cat("q is",qitrs,"\n") }
                        
						qitrs <- - Z[i,r] / F0 * (qitrs * (F(theta0[s] + theta/qitrs) - F0) - A/(1-A)*theta)
                        if (is.nan(qitrs)) {
                            stop("qitrs is nan \n")
                        }
    					if (qitrs <= 0) { 
                            qitrs <- EPSILON
                            break
                        }
                        
						k <- k + 1
					} # while
                    q[i,t,r,s] <- qitrs
					#cat("i:",i,"t:",t,"s:",s,"r:",r,"q:",q[i,t,r,s],"\n")
				} # r
				
				# normalize w.r.t. findings (Y)
				q[i,t,,s] <- q[i,t,,s] / sum(q[i,t,,s])
			} # s
		} # t
	} # i 
	
	return(q)
}

UpdateVariationalParams.C <- function(params,Z,confusion,Y,Xo,Xd,visits)
{
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	nSpecies <- dim(Y)[3]
	
	if (class(params) != "list") { stop("parmas must be a list...") }
	
	alphas <- betas <- beta0s <- NULL
	for (s in 1:nSpecies) {
		alphas <- c(alphas,params[[s]]$alpha)
		
		beta0s <- c(beta0s, params[[s]]$beta0)
		
		for (r in 1:nSpecies) {
			betas <- c(betas,params[[s]]$beta[[r]])
		}
	}
    
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tZ  <- t(Z)
	
	# call the C subroutine to compute the derivs of EJLL
	result <- .C("UpdateVariationalParams", 
                 as.integer(nSites), 
                 as.integer(nVisits), 
                 as.integer(nOccCovs), 
                 as.integer(nDetCovs), 
                 as.integer(nSpecies),
                 as.double(alphas), 
                 as.double(betas), 
                 as.double(beta0s), 
			     as.double(tZ), 
                 as.integer(tY), 
                 as.double(tXo), 
                 as.double(tXd), 
                 as.integer(visits),
			     q = double(nSites*nVisits*nSpecies*nSpecies) )
	q <- result[["q"]]
	
	# reshape
	#dim(q) <- c(nSites,nVisits,nSpecies,nSpecies)
    #q <- aperm(q,c(1,2,3,4))
	dim(q) <- c(nSpecies,nSpecies,nVisits,nSites)
	q <- aperm(q,c(4,3,2,1))
	
	if (sum(is.nan(q)) > 0) { 
		#print(q)
		stop("q has nan")
	}
	
	return(q)
}


#
# Expectation-Maximization for MSODVL model
# In the VEM, first randomly initialize the expected site occupancies.
# During learning, do the M-step first and then the E-step
#
EM.MSODVL <- function(confusion,Y,Xo,Xd,visits,regType,lambda,leakProbs,trueZ=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	# leak probability
	if (is.null(leakProbs)) { stop("We do not support non-fixed leak probabilities...\n") } 
	fixedLeak <- TRUE
	
	# first E-step: randomly assign probabilities for Z and model parameters
	confusion <- matrix(1,nSpecies,nSpecies)
	nParams <- nOccCovs*nSpecies+nDetCovs*nSpecies*nSpecies
	params <- array(rnorm(nParams,mean=0,sd=0.01), c(nParams, 1))
    #params <- array(rep(0,nParams), c(nParams, 1))
	paramsList <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs)
	
	# initialize Z
	Z <- matrix(rep(0.1,nSites*nSpecies),nrow=nSites)
	for (s in 1:nSpecies) {
		Z[which(rowSums(Y[,,s]) > 0),s] <- 0.9
	}
	if (is.null(trueZ) != TRUE) { Z <- trueZ }
	
	newEJLL <- -Inf
    cat("Initial params:",params,"\n")
	diffParams <- 1/EPSILON
	iteration <- 1	
	while (diffParams > TOLERANCE && iteration <= MAX_ITERATIONS) {
        #cat("Iteration:",iteration,"\n")
        oldParams <- params
        
		# E-step: first update q and then calcualte expected Z
		q <- UpdateVariationalParams.C(paramsList,Z,confusion,Y,Xo,Xd,visits)
        #qq <- UpdateVariationalParams.R(paramsList,Z,confusion,Y,Xo,Xd,visits)
		#cat("=== q is done ===\n")

		Z <- ComputeExpectedOccs.MSODVL.C(paramsList,q,confusion,Y,Xo,Xd,visits)
        #ZZ <- ComputeExpectedOccs.MSODVL.R(paramsList,q,confusion,Y,Xo,Xd,visits)
		#cat("=== Z is done ===\n")
				
		# M-step: update the model parameters
		outputs <- optim(params,ComputeEJLL.MSODVL.C,ComputeDerivsOfEJLL.MSODVL.C,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs,
						method="BFGS",
						control=list(maxit=BFGS_MAX_ITERS)) # BFGS  CG  trace=6 in control
		params <- outputs$par
		#cat("=== theta is done ===\n")
		
		# udpate params
		diffParams <- sum((params-oldParams)^2) / length(params)
		
		# compute LL and EJLL
		paramsList <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs)
		if (iteration %% 5 == 0 && PRINTTRACE) { 
			newEJLL <- -ComputeEJLL.MSODVL.C(paramsList,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs)
			newLL <- 0 # -ComputeLL.MSODVL(paramsList,confusion,Y,Xo,Xd,visits,regType,lambda)
			cat("EM.MSODVL iteration: ", iteration, " EJLL is ", newEJLL, " LL is ", newLL, "params change is ", diffParams, "\n") 
		}
		iteration <- iteration + 1
	} # while
	cat("VEM converges in",iteration,"iterations.\n")
    
	newEJLL <- -ComputeEJLL.MSODVL.C(paramsList,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs)
	return(list(params=paramsList,q=q,Z=Z,EJLL=newEJLL))
}

#
# random restart EM for MSODVL model
#
RandomRestartEM.MSODVL <- function(confusion,Y,Xo,Xd,visits,regType=2,lambda,nRandomRestarts=2,leakProbs,trueZ=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	fixedLeak <- TRUE
	
	# random restart
	bestLL <- bestEJLL <- -Inf
	for (i in 1:nRandomRestarts) {
		cat("---------------------------\n")
		cat("MSODVL RandomRestart Iteration",i,"\n")
		
		done <- FALSE
		while (!done) 
		{
			r <- try(result <- EM.MSODVL(confusion,Y,Xo,Xd,visits,regType,lambda,leakProbs,trueZ))
			done <- !inherits(r, "try-error")
		}
		
		newEJLL <- result$EJLL 
		cat("Final EJLL is", newEJLL, "\n")
		newLL <- result$EJLL
		cat("Final LL is", newLL, "\n")
		
		if (newLL > bestLL) {
			bestResult <- result
			bestLL <- newLL
			bestEJLL <- newEJLL
		}
	} # i
	
	cat("**********************\n")
	cat("Best EJLL is", bestEJLL, "\n")
	print(bestResult$params)
	
	return(list(params=bestResult$params,
				confusion=confusion,
				EJLL=bestResult$EJLL,
				LL=bestResult$LL,
				fixedLeak=fixedLeak
				))
}

#
# predict site occupancy
#
PredictOcc.MSODVL <- function(params,confusion,Y,Xo,Xd,visits)
{    
	nSpecies <- nrow(confusion)   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	# convert argument params into list
	if (class(params) != "list") { stop("parmas must be a list...") }
		
	# initialize Z
	Z <- matrix(rep(0.5,nSites*nSpecies),nrow=nSites)
	for (s in 1:nSpecies) {
		Z[which(rowSums(Y[,,s]) > 0),s] <- 0.8
	}
	
	maxIterations <- 10
	iteration <- 1
	while (iteration <= maxIterations) {
		q <- UpdateVariationalParams.C(params,Z,confusion,Y,Xo,Xd,visits)
		Z <- ComputeExpectedOccs.MSODVL.C(params,q,confusion,Y,Xo,Xd,visits)

		iteration <- iteration + 1
	} # while
	
	return(Z)
}

#
# predict observations
#
PredictDet.MSODVL <- function(params,confusion,Xo,Xd,visits) 
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
        result <- .C("PredictDetNoisyOR", 
                     as.integer(nSites), 
                     as.integer(nVisits), 
                     as.integer(nOccCovs), 
                     as.integer(nDetCovs), 
                     as.integer(nSpecies), 
                     as.double(alphas), 
                     as.double(betas), 
                     as.double(beta0s), 
                     as.integer(confusion[,s]), 
                     as.double(Xo), 
                     as.double(Xd), 
                     as.integer(visits), 
                     as.integer(Zs), 
                     as.integer(nZs),
                     predictedDet = double(nSites*nVisits))
		
		predictedDet[,,s] <- t(array(result[["predictedDet"]],c(nVisits,nSites)))
	} # s
	
	return(predictedDet)
}

#
# predict occ and det
#
Predict.MSODVL <- function(params,confusion,Y,Xo,Xd,visits)
{
	probOcc <- PredictOcc.MSODVL(params,confusion,Y,Xo,Xd,visits)
	probDet <- PredictDet.MSODVL(params,confusion,Xo,Xd,visits)
	
	return(list(probOcc=probOcc,probDet=probDet))
}

#####################
# R and C functions #
#####################

#
# compute loglikelihood of species s at site i for MSODVL
# params has to be a list
#
ComputeiLL.MSODVL.R <- function(params,qi,confusion,Yi,Xoi,Xdi,visit) 
{
	nSpecies <- dim(confusion)[1]
	
	if (class(params) != "list") { stop("Params has to be a list.") }
	
	theta0 <- array(0,nSpecies)
	for (s in 1:nSpecies) {
		theta0[s] <- -log(1-params[[s]]$beta0)
	}
	
	Z1 <- Z0 <- array(0,nSpecies)
	for (r in 1:nSpecies) {
		Oir <- Logistic(Xoi %*% params[[r]]$alpha)
		
		Z1[r] <- Oir
		Z0[r] <- 1-Oir
		for (t in 1:visit) {
			for (s in 1:nSpecies) {
				ditrs <- Logistic(Xdi[t,] %*% params[[s]]$beta[[r]])
				if (ditrs == 0) { ditrs <- EPSILON }
				if (ditrs == 1) { ditrs <- 1-EPSILON }
				theta <- -log(1 - ditrs) 
				
				if (Yi[t,s] == 1) {
					probDet1 <- exp(qi[t,r,s]*(F(theta0[s]+theta/qi[t,r,s]) - F(theta0[s])) + qi[t,r,s]*F(theta0[s]))
					probDet0 <- exp(qi[t,r,s]*F(theta0[s]))
				} else {
					probDet1 <- exp(- theta)
					probDet0 <- 1
				}
				
				Z1[r] <- Z1[r] * probDet1
				Z0[r] <- Z0[r] * probDet0
			} # s
		} # t
	} # r
	
#	return(list(Z1=Z1,Z0=Z0))
	return(Z1/(Z1+Z0))
}

ComputeiLL.MSODVL.C <- function(params,q,confusion,Yi,Xoi,Xdi,visit) 
{
	nOccCovs <- length(Xoi)
	nDetCovs <- dim(Xdi)[2]
	nSpecies <- dim(confusion)[1]
	
	Zs  <- c(0,1)
	nZs <- nrow(Zs)
	
	if (class(params) != "list") { stop("Params has to be a list.") }
	
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
	result <- .C("ComputeiLLMSODVL", as.integer(nOccCovs), as.integer(nDetCovs), as.integer(nSpecies),
			as.double(alphas), as.double(betas), as.double(beta0s), as.integer(confusion), as.integer(tYi),
			as.double(Xoi), as.double(tXdi), as.integer(visit), as.integer(tZs), as.integer(nZs), 
			Z = double(2*nSpecies))
	
	Z <- result[["Z"]]
	return(list(Z1=Z[1:nSpecies],Z0=Z[(nSpecies+1):(2*nSpecies)]))
}


#
# compute expected joint log-likelihood for MSODVL
# the params has to be array
#
ComputeEJLL.MSODVL.R <- function(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	lambdaO <- lambda$O
	lambdaD <- lambda$D
	
	if (class(params) != "list") {
		params <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
	# compute site occupancy probability for each species
	O <- array(0, c(nSites, nSpecies))
	theta0 <- array(0,nSpecies)
	theta <- array(0,c(nSites,nVisits,nSpecies,nSpecies))
	for (s in 1:nSpecies) {
		O[,s] <- Logistic(Xo %*% params[[s]]$alpha)
		
		theta0[s] <- -log(1-params[[s]]$beta0)
		
		for (r in 1:nSpecies) {
			for (t in 1:nVisits) {
				theta[,t,r,s] <- -log(1-Logistic(Xd[,t,] %*% params[[s]]$beta[[r]]))
			} # t
		} # r
	} # s
	O[O == 0] <- EPSILON
	O[O == 1] <- 1-EPSILON
	theta[is.infinite(theta)] <- 1/EPSILON
	
	EJLL <- 0
	for (i in 1:nSites) {     
		Xdi <- matrix(Xd[i,,],nrow=nVisits)
		visit <- visits[i]
		iEJLL <- 0
		
		# term on occupancy
		for (r in 1:nSpecies) {
			iEJLL <- iEJLL + Z[i,r] * log(O[i,r]/(1-O[i,r])) + log(1-O[i,r])
		} # r
		
		# term on detection
		for (t in 1:visit) {
			for (s in 1:nSpecies) {
				iEJLL <- iEJLL + (- theta0[s] - Z[i,]%*%theta[i,t,,s]) * (1-Y[i,t,s])
				
				for (r in 1:nSpecies) {
					iEJLL <- iEJLL + Z[i,r] * Y[i,t,s] * q[i,t,r,s] * (F(theta0[s] + theta[i,t,r,s]/q[i,t,r,s]) - F(theta0[s])) + Y[i,t,s] * q[i,t,r,s] * F(theta0[s])
				} # r
			} # s
		} # t
		iEJLL <- as.numeric(iEJLL)
		
		if (iEJLL == Inf || iEJLL == -Inf) { stop("Site ",i," iEJLL is Inf.\n") }
		if (is.na(iEJLL)) { stop("Site ",i," iEJLL is NA.\n") }
		if (is.nan(iEJLL)) { stop("Site ",i," iEJLL is NaN.\n") }
		
		EJLL <- EJLL + iEJLL
	} # i
		
	# regularization
	if (regType > 0) {
		for (r in 1:nSpecies) {
			EJLL <- EJLL - 0.5*lambdaO*sum(params[[r]]$alpha[2:nOccCovs]^2) 
			
			for (s in 1:nSpecies) {
				EJLL <- EJLL - 0.5*lambdaD*sum(params[[s]]$beta[[r]][2:nDetCovs]^2)
			} # s
		} # r
	}
	
	# penalty
	if (PENALTY == 1) {
		Xdd <- NULL
		for (i in 1:nSites) {
			Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
		}
		for (r in 1:nSpecies) {
			avgTrueDetProb <- mean(Logistic(Xdd %*% params[[r]]$beta[[r]]))
			for (s in 1:nSpecies) {
				if (r == s) { next }
				avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[r]]))
				if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
					EJLL <- EJLL - LAMBDAI * (avgFalseDetProb - avgTrueDetProb)
				}
			}
		}
	}
	
	negEJLL <- -EJLL
	if (EJLL == Inf || EJLL == -Inf) { stop("EJLL is Inf.\n") }
	if (is.na(EJLL)) { stop("EJLL is NA...\n") }
	if (is.nan(EJLL)) { stop("EJLL is NaN...\n") } 
	
	return(negEJLL)
}

ComputeEJLL.MSODVL.C <- function(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]  
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	
	lambdaO <- lambda$O
	lambdaD <- lambda$D
	
	if (class(params) != "list") {
		params <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
	alphas <- betas <- beta0s <- NULL
	for (s in 1:nSpecies) {
		alphas <- c(alphas,params[[s]]$alpha)
		
		beta0s <- c(beta0s, params[[s]]$beta0)
		
		for (r in 1:nSpecies) {
			betas <- c(betas,params[[s]]$beta[[r]])
		}
	}
    
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tZ  <- t(Z)
	tq  <- aperm(q,c(4,3,2,1))
	
	# call the C subroutine to compute expected joint log-likelihood
	result <- .C("ComputeEJLLMSODVL", 
				as.integer(nSites), 
				as.integer(nVisits), 
				as.integer(nOccCovs), 
				as.integer(nDetCovs), 
				as.integer(nSpecies),
				as.double(alphas), 
				as.double(betas), 
				as.double(beta0s), 
				as.double(tq),
				as.double(tZ), 
				as.integer(tY), 
				as.double(tXo), 
				as.double(tXd), 
				as.integer(visits), 
				EJLL = double(1))
	EJLL <- result[["EJLL"]]
	
	# regularization
	if (regType > 0) {
		for (r in 1:nSpecies) {
			EJLL <- EJLL - 0.5*lambdaO*sum(params[[r]]$alpha[2:nOccCovs]^2) 
			
			for (s in 1:nSpecies) {
				EJLL <- EJLL - 0.5*lambdaD*sum(params[[s]]$beta[[r]][2:nDetCovs]^2)
			} # s
		} # r
	}
	
	# penalty
	if (PENALTY == 1) {
		Xdd <- NULL
		for (i in 1:nSites) {
			Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
		}
		for (r in 1:nSpecies) {
			avgTrueDetProb <- mean(Logistic(Xdd %*% params[[r]]$beta[[r]]))
			for (s in 1:nSpecies) {
				if (r == s) { next }
				avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[r]]))
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
# compute derivatives of parameter w.r.t. EJLL for MSODVL
# the params has to be an array
#
ComputeDerivsOfEJLL.MSODVL.R <- function(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]   # number of species
	nSites   <- dim(Xd)[1]  # number of sites
	nVisits  <- dim(Xd)[2]  # number of visits
	nOccCovs <- dim(Xo)[2]  # number of occupancy covs
	nDetCovs <- dim(Xd)[3]  # number of detection covs
	
	lambdaO <- lambda$O
	lambdaD <- lambda$D
	
	if (class(params) != "list") {
		params <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
	# compute site occupancy probability for each species
	O <- array(0, c(nSites, nSpecies))
	theta0 <- array(0,nSpecies)
	theta <- array(0,c(nSites,nVisits,nSpecies,nSpecies))
	for (s in 1:nSpecies) {
		O[,s] <- Logistic(Xo %*% params[[s]]$alpha)
		
		theta0[s] <- -log(1-params[[s]]$beta0)
		
		for (r in 1:nSpecies) {
			for (t in 1:nVisits) {
				theta[,t,r,s] <- -log(1-Logistic(Xd[,t,] %*% params[[s]]$beta[[r]]))
			} # t
		} # r
	} # s
	O[O == 0] <- EPSILON
	O[O == 1] <- 1-EPSILON
	theta[is.infinite(theta)] <- 1/EPSILON
	
	dQdalpha <- list()
	dQdbeta <- list()
	for (s in 1:nSpecies) {
		dQdalpha[[s]] <- rep(0,nOccCovs)
		
		dQdbeta[[s]] <- list()
		for (r in 1:nSpecies) {
			dQdbeta[[s]][[r]] <- rep(0,nDetCovs)
		} # ss
	} # s
	
	# compute derivatives
	for (s in 1:nSpecies) {
		# derivs of alpha
		for (i in 1:nSites) {
			dQdalpha[[s]] <- dQdalpha[[s]] + (Z[i,s]-O[i,s])*Xo[i,]
		} # i
		
		# derivs of beta
		for (r in 1:nSpecies) {
			for (i in 1:nSites) {
				for (t in 1:visits[i]) {
					d <- Logistic(Xd[i,t,] %*% params[[s]]$beta[[r]])
					A <- exp(- theta0[s] - theta[i,t,r,s] / q[i,t,r,s])
					dQdbeta[[s]][[r]] <- dQdbeta[[s]][[r]] + Z[i,r] * (Y[i,t,s]/(1-A) -1) * d *Xd[i,t,]
				} # t
			} # i
		} # r
	} # s
	
	if (regType > 0) {
		for (s in 1:nSpecies) {
			dQdalpha[[s]] <- dQdalpha[[s]] - lambdaO * c(0, params[[s]]$alpha[2:nOccCovs])
			
			for (r in 1:nSpecies) {
				dQdbeta[[s]][[r]]  <- dQdbeta[[s]][[r]] - lambdaD * c(0, params[[s]]$beta[[r]][2:nDetCovs])
			} # ss
		} # s
	}
	
	# extract the derivs
	derivs <- NULL
	for (s in 1:nSpecies) {
		# alpha
		derivs <- c(derivs, dQdalpha[[s]])
		
		# beta
		for (r in 1:nSpecies) {    
			derivs <- c(derivs, dQdbeta[[s]][[r]])
		} # ss
	} # s
	
	derivs <- -derivs
	return(derivs)
}

ComputeDerivsOfEJLL.MSODVL.C <- function(params,q,confusion,Z,Y,Xo,Xd,visits,regType,lambda,fixedLeak,leakProbs=NULL) 
{
	nSpecies <- dim(Y)[3]  
	nSites   <- dim(Xd)[1]
	nVisits  <- dim(Xd)[2]
	nOccCovs <- dim(Xo)[2]
	nDetCovs <- dim(Xd)[3]
	
	lambdaO <- lambda$O
	lambdaD <- lambda$D
	
	if (class(params) != "list") {
		params <- ParamsArray2List.MSODVL(params,confusion,nSpecies,nOccCovs,nDetCovs,fixedLeak,leakProbs) 
	}
	
	alphas <- betas <- beta0s <- NULL
	for (s in 1:nSpecies) {
		alphas <- c(alphas,params[[s]]$alpha)
		
		beta0s <- c(beta0s, params[[s]]$beta0)
		
		for (r in 1:nSpecies) {
			betas <- c(betas,params[[s]]$beta[[r]])
		}
	}
	
	# transpose the matrices
	tY  <- aperm(Y,c(3,2,1))
	tXo <- t(Xo)
	tXd <- aperm(Xd,c(3,2,1))
	tZ  <- t(Z)
	tq  <- aperm(q,c(4,3,2,1))
	
	# call the C subroutine to compute expected joint log-likelihood
	result <- .C("ComputeDerivsOfEJLLMSODVL", 
				as.integer(nSites), 
				as.integer(nVisits), 
				as.integer(nOccCovs), 
				as.integer(nDetCovs), 
				as.integer(nSpecies),
				as.double(alphas), 
				as.double(betas), 
				as.double(beta0s), 
				as.double(tq),
				as.double(tZ), 
				as.integer(tY), 
				as.double(tXo), 
				as.double(tXd), 
				as.integer(visits), 
				derivsOfEJLL = double(nSpecies*nOccCovs+nSpecies*nSpecies*nDetCovs) )
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
			derivs[(count):(count+nOccCovs-1)] <- derivs[(count):(count+nOccCovs-1)] - lambdaO * c(0, params[[s]]$alpha[2:nOccCovs])
			count <- count + nOccCovs
			
			# find the number of violations on detection probability
			numViolations <- 0
			avgTrueDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[s]]))
			for (r in 1:nSpecies) {
				avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[r]]))
				if (avgFalseDetProb > avgTrueDetProb + 1e-10) {
					numViolations <- numViolations + 1
				}
			}
			if (PENALTYALL) {
				numViolations <- nSpecies-1
			}
			
			for (r in 1:nSpecies) {
				if (PENALTY == 1) {
					# penalty
					avgFalseDetProb <- mean(Logistic(Xdd %*% params[[s]]$beta[[r]]))
					falseDetProb <- Logistic(Xdd %*% params[[s]]$beta[[r]])
					multiplier <- falseDetProb * (1-falseDetProb)
					if (r != s) {
						if (avgTrueDetProb < avgFalseDetProb || PENALTYALL) {
							derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] - 
									LAMBDAI * rowMeans(sapply(1:nrow(Xdd), function(i,a,b) a[i]*b[i,], a = multiplier, b = Xdd))
						}
					} else {
						derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] + 
								LAMBDAI * rowMeans(sapply(1:nrow(Xdd), function(i,a,b) a[i]*b[i,], a = multiplier, b = Xdd)) * numViolations
					}
				}
				
				derivs[(count):(count+nDetCovs-1)] <- derivs[(count):(count+nDetCovs-1)] - lambdaD * c(0, params[[s]]$beta[[r]][2:nDetCovs])
				count <- count + nDetCovs
			} # r
		} # s
	}
	
	if (DEBUG) { print(derivs) }
	
	if (sum(is.infinite(derivs)) > 0) { print(derivs); stop("derivs contains Inf.\n") }
	if (sum(is.na(derivs)) > 0)       { print(derivs); stop("derivs contains NA...\n") }
	if (sum(is.nan(derivs)) > 0)      { print(derivs); stop("derivs contains NaN...\n") }  
	
	derivs <- -derivs
	return(derivs)
}
