# Evaluate predicting observations (Y) using different methods on eBird data
#
# Author: Jun Yu
# Version: Oct 2013
# Usage: Rscript MSOD_eBird_eval.R Hawk OD 
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")
source("MSOD.R")
source("MSODL1.R")
source("MSOD_metric.R")

# commandline args
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]  # Hawk Woodpecker Finch
model <- args[2] # OD ODLP MSOD
if (! (dataset %in% c("Hawk","Woodpecker","Finch"))) { stop("invalid dataset...\n") }
if (! (model %in% c("OD","ODLP","MSOD"))) { stop("invalid model...\n") }
cat("===",model,"on",dataset,"===\n")

acc <- auc <- f1 <- NULL

exps <- 11:15
nExps <- length(exps)
for (i in 1:nExps) {
    # load result
    resultFile <- paste("results/eBird/",dataset,"_",model,"_",i,".RData",sep="")
    load(resultFile)
	
	#print(params)
    
    # evaluate
    metric <- evaluateDet(teY,teYHat,teVisits,"True")                    
    
    acc <- rbind(acc,metric$acc)
    auc <- rbind(auc,metric$auc)
    f1  <- rbind(f1,metric$f1)
} # i

for (i in 1:2) {
    cat("- Species",i,"-\n")
#    cat("acc:",mean(sort(acc[,i],decreasing=TRUE)[1:10]),"+",stderr(acc[,i]),"\n")
#    cat("auc:",mean(sort(auc[,i],decreasing=TRUE)[1:10]),"+",stderr(auc[,i]),"\n")
#    cat("f1: ",mean(sort(f1[,i],decreasing=TRUE)[1:10]), "+",stderr(f1[,i]), "\n")
	cat("acc:",mean(acc[,i]),"+",stderr(acc[,i]),"\n")
	cat("acc:",mean(auc[,i]),"+",stderr(auc[,i]),"\n")
	cat("acc:",mean(f1[,i]),"+",stderr(f1[,i]),"\n")
}
#print(acc)
#print(auc)
#print(f1)






