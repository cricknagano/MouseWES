## PER_PATIENT_SCRIPT_PRE_PROCESSING

library(xlsx)

source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/returnpathstofastqfilesformouse.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/Screatebwamemcommandformousefastqfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/ScreateRGgroupsforsamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/ScreateCleanSamcommandforsamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createSortcommandforbamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createBAIfileforbamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createFLAGSTATcommandforfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createMERGEcommandformousebamfiles.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createBAIfileforbamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createDEDUPcommandsforfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createTARGETrealignmentintervallist.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createINDELREALIGNMENTcommands.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createPICARDQCcommandforbamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createFASTQCcommandforbamfile.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createCOVERAGEDEPTHcommands.R')
source('/camp/lab/swantonc/working/bakkerb/exome_seq_mouse_pipeline/auxiliary_functions/createMPILEUPfilesforbam.R')

# set paths (for output etc)
pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/DN19306all.xlsx"
path          <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/fastq/"

#pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"
#path <- "/camp/stp/sequencing/inputs/instruments/fastq/200720_K00102_0489_AHH73WBBXY/fastq/DN19306/"




mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)
mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
#paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/S0276129_Paddedliftovermm10.bed"
paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"

## specify NGSIDs
NGSIDs <- mousedata$Genomics.ID


for(i in 1:length(NGSIDs)){
        
        
        ###################################        LOCATE FASTQ FILES/CREATE ANALYSIS FOLDER        ######################################### 
        

        fastqpathsformouse <- returnpathstofastqfilesformouse(path, NGSIDs[i])
        print(fastqpathsformouse)
        listoffailedmice <- NULL
        mouse.name <- mousedata$Mouse.Name[mousedata$Genomics.ID%in%NGSIDs[i]]
        mouse.name <- unique(mouse.name)
        if(length(fastqpathsformouse)<=1){
                print(paste("no fastq files for: ",  NGSIDs[i], sep=""))
                listoffailedmice <- c(listoffailedmice, NGSIDs[i])
		next()
        }else{
                Dir <- pathtosamfiles
                mouseDir <- paste(Dir, mouse.name, sep="")
                if(dir.exists(mouseDir)){
                        print(paste("Directory for mouse ", NGSIDs[i], " exists"))
                        if(!dir.exists(paste(mouseDir, "/errorlogs/", sep =""))){dir.create(paste(mouseDir, "/errorlogs/", sep =""))} 
                        if(!dir.exists(paste(mouseDir, "/flagstats/", sep =""))){dir.create(paste(mouseDir, "/flagstats/", sep ="")) }
                        if(!dir.exists(paste(mouseDir, "/dedupmetrics/", sep =""))){dir.create(paste(mouseDir, "/dedupmetrics/", sep ="")) }
                        if(!dir.exists(paste(mouseDir, "/PICARDQC/", sep =""))){dir.create(paste(mouseDir, "/PICARDQC/", sep ="")) }
                        if(!dir.exists(paste(mouseDir, "/FASTQC/", sep =""))){dir.create(paste(mouseDir, "/FASTQC/", sep ="")) }
                        if(!dir.exists(paste(mouseDir, "/GATKcoverage/", sep =""))){dir.create(paste(mouseDir, "/GATKcoverage/", sep =""))}
                }else{
                        dir.create(mouseDir, mode = "775")
                        dir.create(paste(mouseDir, "/errorlogs/", sep =""), mode = "775") 
                        dir.create(paste(mouseDir, "/flagstats/", sep =""), mode = "775") 
                        dir.create(paste(mouseDir, "/dedupmetrics/", sep =""), mode = "775") 
                        dir.create(paste(mouseDir, "/PICARDQC/", sep =""), mode = "775") 
                        dir.create(paste(mouseDir, "/ASTQC/", sep =""), mode = "775") 
                        dir.create(paste(mouseDir, "/GATKcoverage/", sep =""), mode = "775")
                        print(paste("Directory for mouse ", NGSIDs[i], " created"), mode = "775")
                }
        }
        
        
        ###################################                 BWA MEM                                 #########################################
                
        bwamemcommands <- Screatebwamemcommandformousefastqfile(fastqpathsformouse, NGSIDs[i], mouseDir)

        
        
        ###################################                 ADD RG SCRIPT                           #########################################
        ##get the names of the sam files that the BWA mem command would have created
        count<- 1
        outputfilevector <- NULL
        while(count < length(fastqpathsformouse)){
                outputfilename <- gsub(paste(".*", NGSIDs[i], sep=""), "", x=fastqpathsformouse[count])
                outputfilename  <- gsub("_R.*", "", x=outputfilename)
                outputfile <- paste(pathtosamfiles, "/", mouse.name , "/", NGSIDs[i], outputfilename,"_", "PEs.sam", sep="")
                outputfilevector <- c(outputfilevector, outputfile)
                count <- count + 2
        }
        samfilespaths <- outputfilevector
        
        
        RGcommands <- ScreateaddRGcommandforsamfile(NGSIDs[i], samfilespaths, mousedatapath)
       
        
        ###################################                 CLEAN SAM FILES                         #########################################

        PEsRGsamfiles <- gsub("PEs.sam", "PEsRG.sam",samfilespaths )

        cleansamcommands <- ScreateCleanSamcommandforsamfile(NGSIDs[i], PEsRGsamfiles)
        

        
        ###################################                 SORTBAM  FILES                          #########################################
        
        PERGbamfiles <- gsub(".sam", ".bam", PEsRGsamfiles)
        
        sortcommands <- createSortcommandforfile(NGSIDs[i], PERGbamfiles)
        
        
        
        ###################################                 INDEX SORTED  FILES                     #########################################
        
        sortedbamfiles <- gsub(".bam", "sorted.bam", PERGbamfiles)
        
        BAIsortedcommands <- createBAIfileforbamfile(NGSIDs[i], sortedbamfiles)
        BAIsortedcommands <-  gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", BAIsortedcommands)
        
        ###################################                 FLAGSTAT SORTED  FILES                  #########################################
        
        flagstatsortedcommands <- createFLAGSTATcommandforfile(NGSIDs[i], sortedbamfiles)
        flagstatsortedcommands <-  gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", flagstatsortedcommands)
        
        ###################################                 MERGE MULTIPLE BAMS                     #########################################

        if (length(sortedbamfiles)==1){
                MERGEcommands <- NULL
        }else{
                MERGEcommands <- createMERGEcommandformousebamfiles(NGSIDs[i], sortedbamfiles)
                MERGEcommands <-  gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", MERGEcommands)
        }
        ###################################                 SORT MERGED FILES                       #########################################
        
        gsublocation <- paste(NGSIDs[i], "_", sep="")
        mergedbamfile <- paste(gsub(paste(gsublocation, ".*", sep=""), "", sortedbamfiles[1]),  NGSIDs[i], "_merged.bam", sep="")
        
        sortmergedcommands <- createSortcommandforfile(NGSIDs[i], mergedbamfile)
        
        ###################################                 INDEX SORTED MERGED BAM                 #########################################
        
        sortedmergedfile <- gsub(".bam", "sorted.bam", mergedbamfile)
        indexmergedsortedfilecommands <- createBAIfileforbamfile(NGSIDs[i], sortedmergedfile)
        
        ###################################                 FLAGSTAT SORTED MERGE BAM FILE          #########################################
        
        flagstatmergedcommands <- createFLAGSTATcommandforfile(NGSIDs[i], sortedmergedfile)
        
        ###################################                 DE_DUP MERGEDSORTED BAM                 #########################################
        
        DEDUPcommands <- createDEDUPcommandsforgivensamfile(NGSIDs[i], sortedmergedfile)
        
        ###################################                 INDEX DEDUP BAM FILE                    #########################################
        
        dedupfile <- gsub(paste(NGSIDs[i], "_", sep=""), paste(NGSIDs[i], "_dedup", sep=""), sortedmergedfile)
        indexdedupcommands <- createBAIfileforbamfile(NGSIDs[i], dedupfile)
        
        ###################################                 FLAGSTAT DEDUP BAMS                     #########################################
        
        flagstatdedupcommands <- createFLAGSTATcommandforfile(NGSIDs[i], dedupfile)
        
        ###################################                 INDELREALIGNER TARGET CREATION          #########################################
        
        targetintervalcommand <- createTARGETrealignmentintervallist(NGSIDs[i], dedupfile, mouserefpath)
        
        ###################################                 INDELREALIGNER REALIGNMENT              #########################################
        
        indelrealignmentcommands <- createINDELREALIGNMENTcommands(NGSIDs[i], dedupfile, mouserefpath)
        
        ###################################                 INDEX PROCESSED FILES                   #########################################
        
        processedfile <- paste( gsub(paste(NGSIDs[i], "_", ".*", sep=""), "", dedupfile), NGSIDs[i], "_processed.bam",  sep="")
        indexprocessedfiles <- createBAIfileforbamfile(NGSIDs[i], processedfile)
        
        ###################################                 CREATE PICARDQC/FASTQC/GATKCOVERAGE     #########################################
        
        PICARDQCcommands <- createPICARDQCcommandforbamfile(NGSIDs[i], processedfile, mouserefpath)
        FASTQCcommands <- createFASTQCcommandforbamfile(NGSIDs[i], processedfile, mouserefpath)
        coveragedepth <- createCOVERAGEDEPTHcommands(NGSIDs[i], processedfile, mouserefpath)
        
        ###################################                 CREATE MPILEUP FILE                     #########################################
        
        mpileupcommands <- createMPILEUPfilesforbam(NGSIDs[i], processedfile, mouserefpath, paddedexonbedfile)
        
        ###################################                 CREATE SHELL SCRIPT                     #########################################
        
        shellDir <- paste(pathtosamfiles, "/",mouse.name, "/", NGSIDs[i], "_pre-processingcommands.sh", sep="")
        
        numberoftasks <- length(bwamemcommands) + length(RGcommands) + length(cleansamcommands) + 
                length(sortcommands) + length(BAIsortedcommands) + length(flagstatsortedcommands) + 
                length(MERGEcommands) + length(sortmergedcommands) + length(indexmergedsortedfilecommands) + 
                length(flagstatmergedcommands) + length(DEDUPcommands) + length(indexdedupcommands) + 
                length(flagstatdedupcommands) + length(targetintervalcommand) + length(indelrealignmentcommands) +
                length(indexprocessedfiles) + length(PICARDQCcommands) + length(FASTQCcommands) +
                length(coveragedepth) + length(mpileupcommands)
        numberofcpus <- 8
        erroroutput <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste(pathtosamfiles, "/", mouse.name, "/errorlogs/", NGSIDs[i], "_", "single_pre-processingcommands.txt", sep=""))
        
        cat("#!/bin/bash", file = shellDir, sep='\n\n', append = FALSE, fill = FALSE)
        
        cat(paste("#SBATCH --time=3-00:00:0"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH -e ", erroroutput, sep=""), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        ##cat(paste("#SBATCH --exclusive", sep=""), file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH --partition=cpu"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH --cpus-per-task=8"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        
        #cat("#SBATCH --mem-per-cpu=7168", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        cat("module load BWA/0.7.15-foss-2016b", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load SAMtools/1.3.1-intel-2016b", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load GATK/3.6-Java-1.8.0_92", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load FastQC/0.10.1-Java-1.7.0_80 ", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load  picard/2.1.1-Java-1.8.0_112 ", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load  R ", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        for (x in 1:length(bwamemcommands)){
                        cat(paste("srun -c 8 -n 1 -J Bmem", x, " time ", bwamemcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
        }
        
        for (x in 1:length(RGcommands)){
                        cat(paste("srun -c 8 -n 1 -J RG", x, " time ", RGcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
        }
        
        for (x in 1:length(cleansamcommands)){
                        cat(paste("srun -c 8 -n 1 -J clean", x, " time ", cleansamcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
        }
        
        for (x in 1:length(sortcommands)){
                        cat(paste("srun -c 8 -n 1 -J sort", x, " time ", sortcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
               
        }
        
        for (x in 1:length(BAIsortedcommands)){
                        cat(paste("srun -c 8 -n 1 -J bai", x, " time ", BAIsortedcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
        }
        
        for (x in 1:length(flagstatsortedcommands)){
                        cat(paste("srun -c 8 -n 1 -J FLGstt", x, " time ", flagstatsortedcommands[x], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                        
        }
        
        if (length(MERGEcommands)>=1){
                cat(paste("srun -c 8 -n 1 -J merge time ", MERGEcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                
        }else{
                gsublocation <- paste(NGSIDs[i], "_", sep="")
                copyfilepath <- paste(gsub(paste(gsublocation, ".*", sep=""), "", sortedbamfiles),  NGSIDs[i], "_merged.bam", sep="")
                cat(gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste("srun -c8 -n 1 -J merge time cp ", sortedbamfiles, " ", copyfilepath, sep="")), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
                
        }
        
        cat(paste("srun -c 8 -n 1 -J sortmerged time ", sortmergedcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J BAImerg time ", indexmergedsortedfilecommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J FLGsttM time ", flagstatmergedcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J dedeup time ", DEDUPcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J BAIdedup time ", indexdedupcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J FLGsttD time ", flagstatdedupcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J interval time ", targetintervalcommand, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J realign time ", indelrealignmentcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J BAIproc time ", indexprocessedfiles, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J PicQC time ", PICARDQCcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J FastQC time ", FASTQCcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J coverage time ", coveragedepth, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J MPIL time ", mpileupcommands, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        
        cat("wait", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)

 }
 
