# Split eBird data into training and test
#
# Author: Jun Yu
# Version: Dec 2012
##################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/MSOD")
setwd("/nfs/guille/tgd/wonglab/yuju/MSOD")

# load data
dataFile <- "./data/eBird_California_2010_msod.RData"
load(dataFile)

Xo     <- siteData$occCovs
Xd     <- siteData$detCovs
Y      <- siteData$detHists
visits <- siteData$visits
nSites <- nrow(Xo)

# split into train and test 
datasets <- list()
for (i in 1:50) {
    allIdx <- 1:nSites
    trIdx <- sample(allIdx,round(nSites/2))
    teIdx <- allIdx[! (allIdx %in% trIdx)]
    
    trXo <- Xo[trIdx,]
    trXd <- Xd[trIdx,,]
    trY  <- Y[trIdx,,]
    trVisits <- visits[trIdx]
    nTrSites <- length(trIdx)
    
    teXo <- Xo[teIdx,]
    teXd <- Xd[teIdx,,]
    teY  <- Y[teIdx,,]
    teVisits <- visits[teIdx]
    nTeSites <- length(teIdx)
 
   datasets[[i]] <- list(trXo=trXo,trXd=trXd,trY=trY,trVisits=trVisits,nTrSites=nTrSites,
                         teXo=teXo,teXd=teXd,teY=teY,teVisits=teVisits,nTeSites=nTeSites)
}

save(datasets,Xo,Xd,Y,visits,nSites,speciesNames,file="./data/eBird_2010_msod.RData")
