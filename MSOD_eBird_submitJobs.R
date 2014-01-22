#! /usr/bin/env Rscript

datasets <- c("Hawk","Woodpecker","Finch")
indices  <- 11:15
models   <- c("OD","ODLP")#,"MSOD")

for (dataset in datasets) 
{
    for (index in indices) 
    {
        for (model in models) 
        {
            script <- paste("#!/bin/csh\n\n", 
                                "#$ -N ", substr(dataset,1,1),".",index,".",model,
                                "\n\n# set working directory on all host to",
                                "\n# directory where the job was started",
                                "\n#$ -cwd",
                                "\n",
                                "\n# send all process STDOUT (fd 2) to this file",
                                "\n#$ -o ./tmp/",dataset,"_",index,"_",model,"_o.txt",
                                "\n",
                                "\n# send all process STDERR (fd) to this file",
                                "\n#$ -e ./tmp/",dataset,"_",index,"_",model,"_e.txt",
                                "\n", 
								"\n#$ -q eecs2",
								"\n", 
                                "\n# Commands \n",
                                "Rscript ./MSOD_eBird_exp.R ",dataset," ",model," ",index,"\n", sep="")
            
            file <- paste("./scripts/",dataset,"_",index,"_",model,".sh", sep="")
            write(script, file=file, append = FALSE)
            system(paste("qsub",file))
            Sys.sleep(1.0)
        }
    }
}
