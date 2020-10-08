library(xlsx)

source('/camp/lab/swantonc/working/albakim/MousePipeline/createSCALPELcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createEXPORTscalpelcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createMUTECTcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createVARSCANcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createPROCESSSOMATICcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createREADCOUNTcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createFPFILTERcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createVCFforvarscancommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createSNPfilecommand.R')


pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"
mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)

paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"
mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
dbsnp <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/mgp.v4.snps.dbSNP_chr.vcf"

## specify mice
mice <- as.character(unique(mousedata$Mouse.Name))


for(i in 1:length(mice)){
  
  mouse.name <- mice[i]
  GERMLINE <- mousedata$Sample.Name[mousedata$sample.type%in%'germline' & mousedata$Mouse.Name%in%mouse.name]
  TUMOURs <- mousedata$Sample.Name[mousedata$sample.type%in%'tumour' & mousedata$Mouse.Name%in%mouse.name]
  
  germlinebam <- paste(pathtosamfiles, "/", mouse.name, "/",  GERMLINE, "_processed.bam", sep="")
  tumourbam <- paste(pathtosamfiles, "/", mouse.name, "/",  TUMOURs, "_processed.bam", sep="")
  
  ############################          CHECK VARSCAN/MUTECT/SCALPEL FOLDERS CREATED             ################################
  VARSCANmain <- paste(pathtosamfiles, "/", mouse.name, "/VARSCAN/", sep="")
  MUTECTmain <- paste(pathtosamfiles, "/", mouse.name, "/MUTECT/", sep="")
  SCALPELmain <- paste(pathtosamfiles, "/", mouse.name, "/SCALPEL/", sep="")
  
  VARSCANdir <- paste(pathtosamfiles, "/", mouse.name, "/VARSCAN/", TUMOURs, "vs", GERMLINE, sep="")
  MUTECTdir <- paste(pathtosamfiles, "/", mouse.name, "/MUTECT/", TUMOURs, "vs", GERMLINE, sep="")
  SCALPELdir <- paste(pathtosamfiles, "/", mouse.name, "/SCALPEL/", TUMOURs, "vs", GERMLINE, sep="")
  
  if(!dir.exists(VARSCANmain)){ dir.create(VARSCANmain)}
  if(!dir.exists(MUTECTmain)){ dir.create(MUTECTmain)}
  if(!dir.exists(SCALPELmain)){ dir.create(SCALPELmain)}
  
  if(any(!dir.exists(VARSCANdir))){ 
    for(j in 1:length(VARSCANdir)){
      dir.create(VARSCANdir[j])
    }
    }
  if(any(!dir.exists(MUTECTdir))){ 
    for(j in 1:length(MUTECTdir)){
      dir.create(MUTECTdir[j])
    }
    }
  if(any(!dir.exists(SCALPELdir))){ 
    for(j in 1:length(SCALPELdir)){
      dir.create(SCALPELdir[j])
    }
  }
  
  #############################          CREATE SCALPEL COMMANDs                                  ################################
  
  scalpelcommand <- createSCALPELcommand(TUMOURs, GERMLINE, tumourbam, germlinebam, mouserefpath, paddedexonbedfile, SCALPELdir)
  scalpelexportcommands <- createEXPORTscalpelcommand(TUMOURs, GERMLINE, tumourbam, germlinebam, mouserefpath, paddedexonbedfile, SCALPELdir)
  
  #############################          CREATE MUTECT COMMANDs                                  ################################
  
  mutectcommand <- createMUTECTcommand(TUMOURs, GERMLINE, tumourbam, germlinebam, mouserefpath, paddedexonbedfile, MUTECTdir, dbsnp)
  
  #############################          CREATE VARSCAN COMMANDs                                  ################################
  
  germline.mpileup <- gsub(".bam", ".mpileup", germlinebam)
  tumour.mpileup <- gsub(".bam", ".mpileup", tumourbam)
  varscancommand <- createVARSCANcommand(TUMOURs, GERMLINE, tumour.mpileup, germline.mpileup, mouserefpath, paddedexonbedfile,VARSCANdir)
  
  ##process somatic
  varscan.output.base.snp <- paste(VARSCANdir, "/", TUMOURs, "vs", GERMLINE, ".snp", sep="")
  varscan.output.base.indel <- paste(VARSCANdir, "/", TUMOURs, "vs", GERMLINE, ".indel", sep="")
  processsomaticsnp <- createPROCESSSOMATICcommand(TUMOURs, varscan.output.base.snp, VARSCANdir)
  processsomaticindel <- createPROCESSSOMATICcommand(TUMOURs, varscan.output.base.indel, VARSCANdir)
  
  ##bamreadcount
  pathtosomaticfilepostPROCSOM <- paste(varscan.output.base.snp, ".Somatic", sep="")
  positionsfile <- paste(pathtosomaticfilepostPROCSOM, ".bamreadcount.input", sep="" )
  createpositionsfilecommand <- paste("grep Somatic", pathtosomaticfilepostPROCSOM ,"| grep -v _ |","awk '{print $1\"\\t\"$2\"\\t\"$2}' >",
                                      positionsfile)
  createpositionsfilecommand <- gsub(pattern = "/Volumes/lab-swantonc/", replacement = "/camp/lab/swantonc/", createpositionsfilecommand)
  bamreadcountcommand <- createREADCOUNTcommand(TUMOURs, GERMLINE, tumourbam, mouserefpath, positionsfile)
  
  ##fpfilter (for snp files)
  FPfilter.input <- paste(pathtosomaticfilepostPROCSOM, ".fpfilter.input", sep="")
  FPfilt.inputcommand <- paste("grep Somatic", pathtosomaticfilepostPROCSOM,"| grep -v _ |","awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}' >",
                               FPfilter.input);
  FPfilt.inputcommand <- gsub(pattern = "/Volumes/lab-swantonc/", replacement = "/camp/lab/swantonc/", FPfilt.inputcommand)
  bamreadoutcountoutput <- gsub(".input", ".output", positionsfile)
  fpfiltercommand0.05 <- createFPFILTERcommand(TUMOURs,FPfilter.input, bamreadoutcountoutput, 0.05)
  fpfiltercommand0.02 <- createFPFILTERcommand(TUMOURs,FPfilter.input, bamreadoutcountoutput, 0.02)
  
  ## create vcf file varscanoutput (for snp)
  FPsnv0.05 <- paste(varscan.output.base.snp, ".Somatic.fpfilter.output.VAF0.05", sep="")
  Somaticsnv <- paste(varscan.output.base.snp, ".Somatic", sep="")
  VCFsnv <- paste(varscan.output.base.snp, ".SomaticnativeandFP.vcf", sep="")
  snv2vcfcommand <- createVCFforvarscancommand(FPsnv0.05, Somaticsnv, VCFsnv, " 0 0 > ")
  SNVbgzipcommand <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste("bgzip ",VCFsnv, sep=''))
  tabixSNVcommand <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste( "tabix -p vcf ", VCFsnv, ".gz", sep=""))
  
  ## create vcf file varscanoutput (for indel)
  HCindel0.05 <- paste(varscan.output.base.indel, ".Somatic.hc", sep="")
  Somaticindel <- paste(varscan.output.base.indel, ".Somatic", sep="")
  VCFindel <- paste(varscan.output.base.indel, ".Somaticnativeandhc.vcf", sep="")
  indel2vcfcommand <- createVCFforvarscancommand(HCindel0.05, Somaticindel, VCFindel, " 0 1 > ")
  INDELbgzipcommand <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/",paste("bgzip ",VCFindel,sep=''))
  tabixINDELcommand <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/",paste( "tabix -p vcf ", VCFindel, ".gz", sep=""))
  
  ## create SNP file from the 'GL'
  germlinesnpcommand <- createSNPfilecommand(GERMLINE, germline.mpileup, VARSCANdir)
  
  #############################          CREATE SHELL SCRIPT                                      ################################
  
  for(j in 1:length(TUMOURs)){
    t <- TUMOURs[j]
    shellDir <- paste(pathtosamfiles, "/",mouse.name, "/", TUMOURs[j], "vs", GERMLINE, "_SNVandINDELscript.sh", sep="")
    erroroutput <- gsub("/Volumes/lab-swantonc/", "/camp/lab/swantonc/", paste(pathtosamfiles, "/", mouse.name, "/errorlogs/", t, "_", "SNVandINDELerrors.txt", sep=""))
    
    cat("#!/bin/bash", file = shellDir, sep='\n\n', append = FALSE, fill = FALSE)
    
    cat(paste("#SBATCH --time=3-00:00:0"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
    cat(paste("#SBATCH -e ", erroroutput, sep=""), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
    ##cat(paste("#SBATCH --exclusive", sep=""), file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
    cat(paste("#SBATCH --partition=cpu"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
    cat(paste("#SBATCH --cpus-per-task=16"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
    cat("module load Perl/5.24.0-foss-2016b", file = shellDir, sep='\n', append = TRUE)
    #cat("module load scalpel/0.5.3", file = shellDir, sep='\n', append = TRUE)
    cat("module load scalpel", file = shellDir, sep='\n', append = TRUE) # loads scalpel 0.5.4 instead, fixed by Chris Hadjigeorgiou
    cat("module load VarScan/2.4.1-Java-1.7.0_80", file = shellDir, sep='\n', append = TRUE)
    cat("module load MuTect/1.1.7-Java-1.7.0_80", file = shellDir, sep='\n', append = TRUE)
    cat("module load bam-readcount/0.7.4-foss-2016b", file = shellDir, sep='\n', append = TRUE)
    cat("module load HTSlib/1.3.1-foss-2016b", file = shellDir, sep='\n', append = TRUE)
    cat("module load Perl-bundle-TracerX/0.3-foss-2016b-Perl-5.24.0", file = shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    #cat("#SBATCH --mem-per-cpu=7168", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
    
    cat(paste("srun -c 8 -n 1 -J varscan", " time ", varscancommand[grepl(t, varscancommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J SProcSom", " time ", processsomaticsnp[grepl(t, processsomaticsnp)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J IProcSom", " time ", processsomaticindel[grepl(t, processsomaticindel)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J InBAMPos", " time ", createpositionsfilecommand[grepl(t, createpositionsfilecommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J OutBAMPos", " time ", bamreadcountcommand[grepl(t, bamreadcountcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J FPinp", " time ", FPfilt.inputcommand[grepl(t, FPfilt.inputcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J FP0.05", " time ", fpfiltercommand0.05[grepl(t, fpfiltercommand0.05)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J FP0.02", " time ", fpfiltercommand0.02[grepl(t, fpfiltercommand0.02)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J SNVVCF", " time ", snv2vcfcommand[grepl(t, snv2vcfcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J INDVCF", " time ", indel2vcfcommand[grepl(t, indel2vcfcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J SNVzip", " time ", SNVbgzipcommand[grepl(t, SNVbgzipcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J SNVtabx", " time ", tabixSNVcommand[grepl(t, varscancommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J INDELzip", " time ", INDELbgzipcommand[grepl(t, INDELbgzipcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J INDELtabx", " time ", tabixINDELcommand[grepl(t, tabixINDELcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 16 -n 1 -J scalp", " time ", scalpelcommand[grepl(t, scalpelcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J export", " time ", scalpelexportcommands[grepl(t, scalpelexportcommands)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    cat(paste("srun -c 8 -n 1 -J mutect", " time ", mutectcommand[grepl(t, mutectcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE)
    
    #cat(paste("srun -c 8 -n 1 -J snps", " time ", germlinesnpcommand[grepl(t, germlinesnpcommand)], sep=""), file=shellDir, sep='\n', append = TRUE)
    #cat(" ", file = shellDir, sep='\n', append = TRUE) # did not refer to any of the tumour IDs
    
    cat(paste("srun -c 8 -n 1 -J snps", " time ", germlinesnpcommand, sep=""), file=shellDir, sep='\n', append = TRUE)
    cat(" ", file = shellDir, sep='\n', append = TRUE) 
    
    cat("wait", file = shellDir, sep='\n', append = TRUE)
    
  }
  
}

