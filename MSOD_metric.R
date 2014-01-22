# MSOD evaluation functions
# 
# Author: Jun Yu
# Version: Oct, 2012
###############################################################################

library("pROC")

# threshold used to conver probabilities into 1/0
THRESHOLD <- 0.5

stderr <- function(x) sqrt(var(x)/length(x))

#
# compute accuracy. Note the predictions must be thresholde into 0/1
#
metricAccuracy <- function(Y,Z) 
{
	if (length(Y) != length(Z)) {
		stop("Length of Y and Z does not match...")
	}
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	return(sum(Y == Z)/(length(Y)))
}

#
# compute F1. Note the predictions must be thresholde into 0/1
#
metricPreRecF1 <- function(Y,Z) 
{
	if (length(Y) != length(Z)) {
		stop("Length of Y and Z does not match...")
	}
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	# for each label, compute precision, recall and F1
	ZandY <- Y & Z
	precision <- sum(ZandY) / sum(Z)
	recall <- sum(ZandY) / sum(Y)
	F1 <- 2*precision*recall / (precision+recall)
	if (is.nan(F1)) {
		F1 <- 0
	}
	
	return(list(precision=precision,recall=recall,F1=F1))
}

#
# compute species-based AUC
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
metricAUC <- function(Y,Z) 
{
	if (length(Y) != length(Z)) {
		stop("Length of Y and Z does not match...")
	}
	
	if (sum(Y) == 0) {
		auc <- 0.5
	} else {
		auc <- auc(Y,Z)
	}
	
	return(auc)
}

#
# evaluate on accuracy, F1 and AUC
# Y is the true labels
# Z is the predicted labels
#
evaluation <- function(Y,Z)
{
	accuracy <- metricAccuracy(Y,Z)
	preRecF1 <- metricPreRecF1(Y,Z)
	auc <- metricAUC(Y,Z)
	
	return(list(accuracy=accuracy, f1=preRecF1$F1, auc=auc))
	# precision=preRecF1$precision, recall=preRecF1$recall,
}

#
# print evaluation of occ
#
evaluateOcc <- function(Y,Z,model="MSOD",printTrace=FALSE)
{
	nSpecies <- ncol(Y)
	acc <- array(0,nSpecies)
    auc <- array(0,nSpecies)
    f1  <- array(0,nSpecies)
    
	if (printTrace) { cat("Evaluation Occ of ",model,":\n") }
	for (s in 1:nSpecies) {
		occMetrics <- evaluation(Y[,s], Z[,s]) 
        acc[s] <- occMetrics$accuracy
        auc[s] <- occMetrics$auc
        f1[s]  <- occMetrics$f1
        
		if (printTrace) { cat("Occ of species",s,": Acc",occMetrics$accuracy," AUC",occMetrics$auc, " F1",occMetrics$f1,"\n") }
	} # s
    
    return(list(acc=acc,auc=auc,f1=f1))
}

#
# print evaluation of det
#
evaluateDet <- function(Y,Z,visits,model="MSOD",printTrace=FALSE)
{
	nSpecies <- dim(Y)[3]
	acc <- array(0,nSpecies)
    auc <- array(0,nSpecies)
    f1  <- array(0,nSpecies)
	
	if (printTrace) { cat("Evaluation Det of ",model,":\n") }
	for (s in 1:nSpecies) {
		detMetrics <- evaluation(convertObs2Array(Y[,,s],visits), convertObs2Array(Z[,,s],visits)) 
        acc[s] <- detMetrics$accuracy
        auc[s] <- detMetrics$auc
        f1[s]  <- detMetrics$f1
        
		if (printTrace) { cat("Occ of species",s,": Acc",detMetrics$accuracy," AUC",detMetrics$auc," F1",detMetrics$f1,"\n") }
	} # s
    
    return(list(acc=acc,auc=auc,f1=f1))
}


#
# convert observation to array
#
convertObs2Array <- function(Y,visits)
{	
	Z <- NULL
	for (i in 1:length(visits)) {
		Z <- c(Z, Y[i,1:visits[i]])
	}
	
	return(Z)
}


#
# convert det occs to matrix
#
convertDetOccs2Matrix <- function(detOccs,visits)
{	
	detOccsMatrix <- NULL
	for (i in 1:length(visits)) {
		detOccsMatrix <- rbind(detOccsMatrix, detOccs[i,1:visits[i],])
	}
	
	return(detOccsMatrix)
}

#
# convert det occs to matrix
#
convertOccOccs2Matrix <- function(occOccs,visits)
{	
	occOccsMatrix <- NULL
	for (i in 1:length(visits)) {
		for (t in 1:visits[i])
		occOccsMatrix <- rbind(occOccsMatrix, occOccs[i,])
	}
	
	return(occOccsMatrix)
}

#
# compute difference of two set of paraemters
#
ComputeParamsDiff <- function(params1,params2,confusion,model) 
{
	nSpecies <- nrow(confusion)
	
	for (s in 1:nSpecies) {
		cat("Species ",s,":\n")
		cat("alpha diff:",sum((params1[[s]]$alpha-params2[[s]]$alpha)^2) / length(params1[[s]]$alpha),"\n")
		
		if (model == "Noisy-OR" || model == "Noisy-OR-Leak" || model == "Noisy-OR-Leak-One") {
			if (model == "Noisy-OR-Leak" || model == "Noisy-OR-Leak-One") {
				cat("beta0 diff on:",sum((params1[[s]]$beta0-params2[[s]]$beta0)^2) / length(params1[[s]]$beta0),"\n")
			}
			
			for (ss in 1:nSpecies) {
				if (confusion[ss,s] == 0) { next }
				cat("beta diff on",ss,":",sum((params1[[s]]$beta[[ss]]-params2[[s]]$beta[[ss]])^2) / length(params1[[s]]$beta[[ss]]),"\n")
			}
		} else if (model == "Additive") {
			cat("beta diff :",sum((params1[[s]]$beta-params2[[s]]$beta)^2) / length(params1[[s]]$beta),"\n")
		} else {
			stop("model is invalid.\n")
		}
	} # s
}

######### multilabel metrics ################




#
# compute Hamming Loss. Note the predictions must be thresholde into 0/1
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
HammingLoss <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	nRow <- nrow(Y)
	nCol <- ncol(Y)
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	return(sum(Y != Z)/(nRow*nCol))
}

HammingLossLabel <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	nRow <- nrow(Y)
	nCol <- ncol(Y)
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	return(colSums(Y != Z)/nRow)
}

#
# compute Ranking Loss.
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
RankingLoss <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	nRow <- nrow(Y)
	siteAUC <- req(0,nRow)
	for (i in 1:nRow) {
		if (sum(Y[i,]) == 0 || mean(Y[i,]) == 1) {
			siteAUC[i] <- 0.5
			next
		}
		siteAUC[i] <- auc(Y[i,],Z[i,])
	}
	
	return(1 - mean(siteAUC))
}


#
# compute Exact match. Note the predictions must be thresholde into 0/1
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
ExactMatch <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	nRow <- nrow(Y)
	nCol <- ncol(Y)
	match <- apply(Y==Z,1,sum)
	
	return(sum(match==nCol) / nRow)
}


#
# compute species-based AUC
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
SpeciesAUC <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	nCol <- ncol(Y)
	speciesAUC <- rep(0,nCol)
	for (j in 1:nCol) {
		if (sum(Y[,j]) == 0) {
			speciesAUC[j] <- 0.5
			next
		}
		speciesAUC[j] <- auc(Y[,j],Z[,j])
	}
	
	return(mean(speciesAUC))
}

SpeciesAUCLabel <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	nCol <- ncol(Y)
	speciesAUC <- rep(0,nCol)
	for (j in 1:nCol) {
		if (sum(Y[,j]) == 0) {
			speciesAUC[j] <- 0.5
			next
		}
		speciesAUC[j] <- auc(Y[,j],Z[,j])
	}
	
	return(speciesAUC)
}



#
# compute Micro F1. Note the predictions must be thresholde into 0/1
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
MicroMeasure <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	# compute precision, recall and F1
	ZandY <- Y & Z
	precision <- sum(ZandY)/sum(Z)
	recall <- sum(ZandY)/sum(Y)
	F1 <- 2 * precision * recall / (precision + recall)
	
	return(list(precision=precision,recall=recall,F1=F1))
}

#
# compute Macro F1. Note the predictions must be thresholde into 0/1
#
# Y is true label matrix N * L where N is # of test instances and L is # of labels
# Z is prediction matrix N * L where N is # of test instances and L is # of labels
#
MacroMeasure <- function(Y,Z) 
{
	if (nrow(Y) != nrow(Z) || ncol(Y) != ncol(Z)) {
		stop("Dim of Y and Z does not match...")
	}
	
	Z[Z < THRESHOLD] <- 0
	Z[Z >= THRESHOLD] <- 1
	
	# for each label, compute precision, recall and F1
	ZandY <- Y & Z
	precision <- apply(ZandY,2,sum) / apply(Z,2,sum)
	recall <- apply(ZandY,2,sum) / apply(Y,2,sum)
	F1 <- 2*precision*recall / (precision+recall)
	
	precision <- precision[!is.nan(precision)]
	recall <- recall[!is.nan(recall)]
	F1 <- F1[!is.nan(F1)]
#	F1[is.nan(F1)] <- 0
	
	return(list(precision=mean(precision),recall=mean(recall),F1=mean(F1)))
}

#
# evaluate different metrics specified by e
#
Evaluate.MSOD <- function(trueLabels, predLabels, metrics=c("HammingLoss")) 
{
	trueLabels <- as.matrix(trueLabels)
	class(trueLabels) <- "numeric"
	predLabels <- as.matrix(predLabels)
	class(trueLabels) <- "numeric"
	
	results <- list()
	
	# HammingLoss
	if (sum(metrics == "HammingLoss") == 1) 
	{
		results[["HammingLoss"]] <- HammingLoss(trueLabels,predLabels)
		results[["HammingLossLabel"]] <- HammingLossLabel(trueLabels,predLabels)
		#cat("Hamming Loss is",results[["HammingLoss"]],"\n")
	}
	
	# ExactMatch
	if (sum(metrics == "ExactMatch") == 1) 
	{
		results[["ExactMatch"]] <- ExactMatch(trueLabels,predLabels)
		#cat("Exact Match is",results[["ExactMatch"]] ,"\n")
	}
	
	# MicroF1
	if (sum(metrics == "MicroF1") == 1) 
	{
		results[["MicroF1"]] <- MicroMeasure(trueLabels,predLabels)$F1
		if (is.nan(results[["MicroF1"]])) {
			results[["MicroF1"]] = 0
		}
		#cat("Micro F1 is",results[["MicroF1"]],"\n")
	}
	
	# MacroF1
	if (sum(metrics == "MacroF1") == 1) 
	{
		results[["MacroF1"]] <- MacroMeasure(trueLabels,predLabels)$F1
		if (is.nan(results[["MacroF1"]])) {
			results[["MacroF1"]] = 0
		}
		#cat("Macro F1 is",results[["MacroF1"]],"\n")
	}
	
	# RankingLoss
	if (sum(metrics == "RankingLoss") == 1) 
	{
		results[["RankingLoss"]] <- RankingLoss(trueLabels,predLabels)
		#cat("Ranking Loss is",results[["RankingLoss"]],"\n")	
	}
	
	# SpeciesAUC
	if (sum(metrics == "SpeciesAUC") == 1) 
	{
		results[["SpeciesAUC"]] <- SpeciesAUC(trueLabels,predLabels)
		results[["SpeciesAUCLabel"]] <- SpeciesAUCLabel(trueLabels,predLabels)
		#cat("Species AUC is",results[["SpeciesAUC"]],"\n")
	}
	
	return(results)
}

stdErr <- function(x) 
{
	sqrt(var(x)/length(x))
}