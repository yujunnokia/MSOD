# Real world experiment 
#
# Author: Jun Yu
# Version: Dec 2012
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/eBird.MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/eBird.MSOD")
#source("MSOD.R")
source("MSODL1.R")
source("MSOD.metric.R")

PATH <- "./data/"

# commandline args
args <- commandArgs(trailingOnly = TRUE)
caseStudyNum <- as.numeric(args[1])

######################
# experiment settings
######################
if (caseStudyNum == 1) {
	caseStudy <- "Grackle" 
} else if (caseStudyNum == 2) {
	caseStudy <- "Raven"
} else if (caseStudyNum == 3) {
	caseStudy <- "Hawk"
} else if (caseStudyNum == 4) {
	caseStudy <- "Woodpecker"
} else if (caseStudyNum == 5) {
	caseStudy <- "Finch"
} else {
	stop("case study number if invalid...\n")
}
#caseStudy <- "Finch" # Grackle  Raven  Hawk  Woodpecker  Finch

if (caseStudy == "Grackle") {
	dataFile <- paste(PATH,"eBird_Louisiana_2010_msod.RData",sep="")
	studySpecies <- c("Quiscalus_quiscula","Quiscalus_major","Quiscalus_mexicanus",
    "Charadrius_vociferus","Cardinalis_cardinalis","Molothrus_aeneus")
    leakProbs <- c(0.04,0.01,0.0,0.19,0.14,0.0)
} else if (caseStudy == "Raven") {
	dataFile <- paste(PATH,"eBird_Texas_2010_msod.RData",sep="")
	studySpecies <- c("Corvus_corax","Corvus_cryptoleucus","Corvus_brachyrhynchos",  "Cathartes_aura")	
    leakProbs <- c(0.0,0.0,0.02,0.11)
} else if (caseStudy == "Hawk") {
	dataFile <- paste(PATH,"eBird_California_2010_Hawk_msod.RData",sep="")
	studySpecies <- c("Accipiter_striatus","Accipiter_cooperii",  "Cathartes_aura")	
    leakProbs <- c(0.0007,0.01275,0.062)
} else if (caseStudy == "Woodpecker") {
	dataFile <- paste(PATH,"eBird_California_2010_msod.RData",sep="")
	studySpecies <- c("Picoides_pubescens","Picoides_villosus",  "Junco_hyemalis")	
    leakProbs <- c(0.01112,0.003854,0.03519)
} else if (caseStudy == "Finch") {
	dataFile <- paste(PATH,"eBird_California_2010_msod.RData",sep="")
	studySpecies <- c("Carpodacus_mexicanus","Carpodacus_purpureus")#,  "Dendroica_coronata","Sayornis_nigricans")	
    leakProbs <- c(0.0847,0.006141)#,0.04371,0.08787)
} else {
	stop("case study is invalid...\n")
}

nRandomRestarts <- 3 # number of random restarts 

# print the experiment settings
cat("=== Experiment Settings ===\n")
cat("caseStudy num is",caseStudy,"\n")
cat("nRandomRestarts is",nRandomRestarts,"\n")

#################
# regularization
#################
# regularization
regType <- 1 # regularization types: 0 for none, 1 for L1, 2 for L2
lambdaO <- 10  # regularization paramters
lambdaD <- 10
lambdaC <-  rep(1,length(studySpecies))
lambda <- list()
lambda$O <- lambdaO
lambda$D <- lambdaD
lambda$C <- lambdaC
cat("=== lambda ===\n")
print(lambda)

#################
# load data
#################
load(dataFile)
speciesIdx <- NULL
for (species in studySpecies) {
	speciesIdx <- c(speciesIdx, which(species == speciesNames))
}

##########################
# get train and test data
##########################
Xo     <- siteData$occCovs
Xd     <- siteData$detCovs
Y      <- siteData$detHists[,,speciesIdx]
#dim(Y) <- c(nrow(Y), ncol(Y), 1)
visits <- siteData$visits

########
# stats
########
if (1 == 0) {
print(dim(Xo))
print(dim(Xd))
print(dim(Y))

YY <- NULL
Xdd <- NULL
for (i in 1:nrow(Xo)) {
    YY <- rbind(YY, Y[i,1:visits[i],])
    Xdd <- rbind(Xdd, Xd[i,1:visits[i],])
}
print(colMeans(Xo))
print(colMeans(Xdd))
print(colMeans(YY))
}

#################
# run experiment
#################
learntParams <- RandomRestartEM.L1.MSOD(Y,Xo,Xd,visits,regType,lambda,nRandomRestarts,leakProbs) 
print(learntParams)
