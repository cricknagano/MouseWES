##PLATYPUS SCRIPT
## Generate SNP calls for germline bam files using platypus, convert MNPs to SNPs, genotype tumour bam files at these SNP positions 

library(xlsx)

source('/camp/lab/swantonc/working/albakim/MousePipeline/createPLATYPUScommand.R')

pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/cell.xlsx"

#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/DN19266cell.xlsx"

#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/DN19306all.xlsx"


#pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"

mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)


mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"
#paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/S0276129_Paddedliftovermm10.bed"
NGSIDs <- as.character(mousedata$Genomics.ID)
dbsnp <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/mgp.v4.snps.dbSNP_chr.vcf"

##define the mice you would like to perform platypus GL SNP calling 
#mice <- as.character(unique(mousedata$Mouse.Name))
GERMLINE <- mousedata$Sample.Name[mousedata$sample.type%in%'germline']
mice <- gsub("Tail", "", GERMLINE)
for (i in 1:length(GERMLINE)){
        germlinebam <- paste(pathtosamfiles, "/", mice[i], "/",  GERMLINE[i], "_processed.bam", sep="")
        ## create Germlinedirectory
        germline.dir <- paste(pathtosamfiles, "/", mice[i], "/GL", sep="")
        platypus.dir <- paste(pathtosamfiles, "/", mice[i], "/PLATYPUS", sep="")
        platypus.copy.dir <- paste(platypus.dir, "/copy_number/", sep="")
        platypus.SNPs.dir <- paste(platypus.dir, "/SNPs/", sep="")
        platypus.indels.dir <- paste(platypus.dir, "/INDELs/", sep="")
        if (!dir.exists(germline.dir)){dir.create(germline.dir)}
        if (!dir.exists(platypus.dir)){dir.create(platypus.dir)}
        if (!dir.exists(platypus.copy.dir)){dir.create(platypus.copy.dir)}
        if (!dir.exists(platypus.SNPs.dir)){dir.create(platypus.SNPs.dir)}
        if (!dir.exists(platypus.indels.dir)){dir.create(platypus.indels.dir)}
        
        genSNP<- 1
        genINDEL <- 0
        minmapqual <- 40
        outputfile <- paste(platypus.SNPs.dir, GERMLINE[i], ".raw.platypus.SNPs.vcf", sep="")
        highqualSNPs <- createPLATYPUScommand(germlinebam, mouserefpath, paddedexonbedfile, genSNP, genINDEL, minmapqual, outputfile)
        
        genSNP<- 0
        genINDEL <- 1
        minmapqual <- 20
        outputfile <- paste(platypus.indels.dir, GERMLINE[i], ".raw.platypus.indels.vcf", sep="")
        highqualINDELs <- createPLATYPUScommand(germlinebam, mouserefpath, paddedexonbedfile, genSNP, genINDEL, minmapqual, outputfile)

        #############################          CREATE SHELL SCRIPT                                      ################################
        
        shellDir <- paste(pathtosamfiles, "/", mice[i], "/", GERMLINE[i], "_Platypusscript.sh", sep="")
        
        numberofcpus <- 8
        erroroutput <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste(pathtosamfiles, "/", mice[i], "/errorlogs/", GERMLINE[i], "_", "platypuserrors.txt", sep=""))
        
        cat("#!/bin/bash", file = shellDir, sep='\n\n', append = FALSE, fill = FALSE)
        
        cat(paste("#SBATCH --time=3-00:00:0"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH -e ", erroroutput, sep=""), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        ##cat(paste("#SBATCH --exclusive", sep=""), file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH --partition=cpu"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(paste("#SBATCH --cpus-per-task=8"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/camp/apps/eb/software/GCCcore/5.4.0/lib64", file = shellDir, sep='\n', append = TRUE)
	cat("module purge", file = shellDir, sep='\n', append = TRUE)
        cat("module load foss/2016b", file = shellDir, sep='\n', append = TRUE)
        cat("module load GLib/2.47.5-foss-2016b", file = shellDir, sep='\n', append = TRUE)
        cat("module load GCCcore/5.4.0", file = shellDir, sep='\n', append = TRUE)

        cat("module load Biopython/1.68-foss-2016b-Python-2.7.12", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load Platypus/0.8.1-foss-2016b-Python-2.7.12", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load VCFtools/0.1.14-foss-2016b-Perl-5.22.1", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat("module load HTSlib/1.3.1-intel-2016b", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        
        cat(paste("srun -c 8 -n 1 -J highQSNPs", " time ", highqualSNPs, " & ",  sep=""), file=shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(" ", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        cat(paste("srun -c 8 -n 1 -J highQIND", " time ", highqualINDELs, " & ", sep=""), file=shellDir, sep='\n', append = TRUE, fill = FALSE)
        cat(" ", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
        
        cat("wait", file = shellDir, sep='\n', append = TRUE, fill = FALSE)
        
}

