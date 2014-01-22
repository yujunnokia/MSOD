#! /usr/bin/env Rscript

datasets <- c("syn","syn-I","syn-NL") #c("syn","syn-I","syn-NL","syn-I-NL")
indices <- 11:30
models  <- c("MSODVL") # c("TRUE","OD","ODLP","MSODTRUE","MSOD") 

for (dataset in datasets) {
    for (index in indices) {
        for (model in models) {
            script <- paste("#!/bin/csh\n\n", 
                            "#$ -N ",dataset,".",model,".",index,
                            "\n\n# set working directory on all host to",
                            "\n# directory where the job was started",
                            "\n#$ -cwd",
                            "\n",
                            "\n# send all process STDOUT to this file",
                            "\n#$ -o ./tmp/",dataset,"_",index,"_",model,"_o.txt",
                            "\n",
                            "\n# send all process STDERR to this file",
                            "\n#$ -e ./tmp/",dataset,"_",index,"_",model,"_e.txt",
                            "\n", 
                            "\n#$ -q eecs1", # eecs
                            "\n", 
                            "\n# Commands \n",
                            "Rscript ./MSOD_syn_exp.R ",dataset," ",index," ",model,"\n", sep="")
        
            file <- paste("./scripts/",dataset,"_",index,"_",model,".sh", sep="")
            write(script, file=file, append = FALSE)
            system(paste("qsub",file))
            Sys.sleep(0.5)
        }
    }
}
