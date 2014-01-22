#! /usr/bin/env Rscript

datasets <- c("syn","syn-I","syn-NL") #c("syn","syn-I","syn-NL","syn-I-NL")
indices <- 1:15
models  <- c("MSODVL") # c("TRUE","OD","ODLP","MSODTRUE","MSOD") 

for (dataset in datasets) {
    for (index in indices) {
        for (model in models) {
            job <- paste("S.",model,".",index, sep="")
        
            system(paste("qdel",job))
        }
    }
}
