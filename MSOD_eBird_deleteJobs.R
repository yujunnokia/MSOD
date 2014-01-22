#! /usr/bin/env Rscript

datasets <- c("Hawk","Woodpecker","Finch")
indices  <- 1:10
models   <- c("OD","ODLP","MSOD")

for (dataset in datasets) 
{
    for (index in indices) 
    {
        for (model in models) 
        {
            job <- paste(substr(dataset,1,1),".",index,".",model, sep="")
            
            system(paste("qdel",job))
        }
    }
}
