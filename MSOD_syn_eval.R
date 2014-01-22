# evaluate different models in predicting occupancy and detection
#
# Author: Jun Yu
# Version: Oct 2013
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSOD.R")
source("MSOD_genData.R")
source("MSOD_metric.R")

# commandline args
args <- commandArgs(trailingOnly = TRUE)
dataType <- args[1]
model    <- args[2]
predType <- args[3] 


if (! (dataType %in% c("syn","syn-I","syn-NL","syn-I-NL"))) { stop("The dataType is invalid...") }
if (! (model %in% c("TRUE","OD","ODLP","MSODTRUE","MSOD","MSODVL"))) { stop("The model is invalid...") }
if (! (predType %in% c("occ","det"))) { stop("The prediction type is invalid...") }
cat("=======================================\n")
cat(" Data type:",dataType,"\n","Model:",model,"\n","predType:",predType,"\n")

exps <- 1:30
#exps <- exps[c(-4)]
nExps <- length(exps)

metrics <- list(acc=array(0,nExps),auc=array(0,nExps),f1=array(0,nExps),strAuc=array(0,nExps))
for (i in 1:nExps) {
    n <- exps[i]
    
    # load dataset (labels)
    if (dataType == "syn") {
        load("./data/MSOD_syndata.RData")
    } else if (dataType == "syn-I") {
        load("./data/MSOD_syndata_nonIndependent.RData")
    } else if (dataType == "syn-NL") {
        load("./data/MSOD_syndata_nonLinear.RData")
    } else {
        load("./data/MSOD_syndata_nonIndependent_nonLinear.RData")
    }
    
    # load prediction
    resultFile <- paste("results/syn/MSOD_",dataType,"_",n,"_",model,"_output.RData",sep="")
    load(resultFile)

    dataset    <- datasets[[n]]
    nTeSites   <- dataset$nTeSites
    teVisits   <- dataset$teVisits
    teOccCovs  <- dataset$teOccCovs
    teDetCovs  <- dataset$teDetCovs
    teDetHists <- dataset$teDetHists
    teTrueOccs <- dataset$teTrueOccs
    
    result <- output$result
    testOccPrediction <- result$occPrediction
    testDetPrediction <- result$detPrediction
    
    # evaluate
    if (predType == "occ") {
        metric <- evaluateOcc(teTrueOccs,testOccPrediction,model)
    } else {
        metric <- evaluateDet(teDetHists,testDetPrediction,teVisits,model)
    }
            
    metrics$acc[i] <- mean(metric$acc)
    metrics$auc[i] <- mean(metric$auc)
    metrics$f1[i]  <- mean(metric$f1)
    
    # compute auc of the learnt structure
    if (model == "MSOD") {
		trueParams <- dataset$trueParams
        metrics$strAuc[i] <- auc(array(trueParams$confusion),array(output$MSODFULLParams$detRates))
    }
}

cat("auc-mean:",mean(metrics$auc),"auc-se:",stderr(metrics$auc),"\n")
cat("acc-mean:",mean(metrics$acc),"acc-se:",stderr(metrics$acc),"\n")
cat("f1-mean:", mean(metrics$f1), "f1-se:", stderr(metrics$f1),"\n")
if (model == "MSOD") {
    cat("strAuc-mean:",mean(metrics$strAuc),"strAuc-sd:",stderr(metrics$strAuc),"\n")
}
print(metrics$auc)
