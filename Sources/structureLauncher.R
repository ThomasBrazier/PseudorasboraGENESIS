############################################################################
#   INRA - ANR project Pseudorasbora
#   STRUCTURE Launcher
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns


#==========================================================
# GENERATE FILES .SH AND QSUB TO LAUNCH STRUCTURE ON A LINUX SERVER
#==========================================================

structureLauncher=function(project="project",path="",niter=c(1,10),Kmin=1,Kmax=10,popsize=727){
  # project="project, the name of the project
  # path="", the path of the final working directory in the server
  # niter=c(1,10), numerotation of multiple successive iterations for the same value of K, necessary for Evanno's method
  # Kmin, minimum number of K to test
  # Kmax, maximum number of K to test
  # popsize=727, population size has to be documented
  
  # Init a launcher object that will contain all commands
  launcher=c()
  # Go to working directory
  launcher=rbind(launcher,
                 paste("cd",path,sep=" "),
                 "",
                 "# Set programs as executables",
                 "chmod +x ./structure",
                 "# Convert data files ot Unix encoding",
                 "# Data files can be recognized by the pattern '00_' e.g. '200_Pparva_structure_native.txt'",
                 "dos2unix *00_*",
                 "",
                 "# Launch batch runs for all values of K and every iterations")

    for (j in Kmin:Kmax) {
      # iteration is set with a random seed between 1000 and 9999, and output of structure is saved in a txt file
      iter=paste("nohup ./structure -K ",j," -D ",sample(1:9999,1)," -o Results/K",j,"run",niter[1]," > Results/outputK",j,"run",niter[1],".txt",sep="")
      for (i in (niter[1]+1):niter[2]) {
        iter=paste(iter,paste(" & nohup ./structure -K ",j," -D ",sample(1:9999,1)," -o Results/K",j,"run",i," > Results/outputK",j,"run",i,".txt",sep=""),sep="")
      }
      launcher=rbind(launcher,iter)
    }
    
    launcher=rbind(launcher,
                   "",
                   "# STRUCTURE Harvester on results to:",
                   "# apply Evanno method (2005) to choose the best K",
                   "# prepare files for CLUMPP",
                   "python structureHarvester.py --dir=./Results/ --out=./Harvester/ --evanno --clumpp",
                   "",
                   "# Get the values of the sampling chain from 'outputK-run-.txt'",
                   "# args are in order Kmin, Kmax and number of replicated runs",
                   "chmod +x ./structure_sampling_chain.sh",
                   paste("./structure_sampling_chain.sh", Kmin, Kmax, niter[2],sep=" "),
                   "",
                   "# CLUMPP on files formated by STRUCTURE Harvester",
                   "# Launch a CLUMPP task on each K (except K=1, no segmentation possible on a single cluster)",
                   "chmod +x ./CLUMPP/CLUMPP",
                   paste("for ((i = ",niter[1],"; i < ",niter[2]+1,"; i++ )); do",sep=""),
                   paste("./CLUMPP/CLUMPP ./CLUMPP/paramfile -i './Harvester/K'$i'.indfile' -p './Harvester/K'$i'.popfile' -o './CLUMPP/K'$i'.outfile' -j './CLUMPP/K'$i'.miscfile' -k $i -c ",popsize," -r ",(niter[2]-niter[1]+1)," &",sep=""),
                   "done",
                   "",
                   "# Distruct to display nice graphs of STRUCTURE, with labels",
                   "# Automatically produce bar plots for each K value",
                   "",
                   "# Output files are copied to the 'data' directory in 'Distruct'",
                   paste("for ((i = ",niter[1],"; i < ",niter[2]+1,"; i++ )); do",sep=""),
                  "cp './CLUMPP/K'$i'.outfile' './Distruct/Data/K'$i'.outfile'",
                   "head -n 18 './Harvester/K'$i'.popfile' >> './Distruct/Data/K'$i'.popfile'",
                  "done",
                   "",
                  "cd Distruct",
                   "chmod +x ./distructLinux1.1",
                   paste("for ((i = ",niter[1],"; i < ",niter[2]+1,"; i++ )); do",sep=""),
                   paste("./distructLinux1.1 -K $i -M 18 -N ",popsize," -p 'Data/K'$i'.popfile' -i 'Data/K'$i'.outfile' -o 'Plots/DistructPlotK'$i'.ps'",sep=""),
                   "done",
                  "",
                  "# Convert .ps as .tiff",
                  "cd Plots/",
                   paste("for ((i = ",niter[1],"; i < ",niter[2]+1,"; i++ )); do",sep=""),
                   "gs -sDEVICE=tiff64nc -r300 -sPAPERSIZE=a4 -dBATCH -dNOPAUSE -sOutputFile='DistructPlotK'$i'.tiff' 'DistructPlotK'$i'.ps'",
                   "done",
                  "cd ..",
                  "cd ..")

      ## Write data file for bash
    if (!dir.exists(paste("Data/STRUCTURE/",project,sep=""))){
      dir.create(paste("Data/STRUCTURE/",project,sep=""))
    }
    file.create(paste("Data/STRUCTURE/",project,"/structure.sh",sep=""),overwrite=TRUE)
    write.table(rbind("#!/bin/sh","# file: structure.sh\n",launcher),paste("Data/STRUCTURE/",project,"/structure.sh",sep=""),
                quote=F,row.names = F,col.names = F)
  }



#==========================================================
# THE END