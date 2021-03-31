##SINGLE_BAMINTERMEDIATECLEANUP_SCRIPT.R

library(xlsx)

source('/camp/lab/swantonc/working/albakim/MousePipeline/returnpathstofastqfilesformouse.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/Screatebwamemcommandformousefastqfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/ScreateRGgroupsforsamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/ScreateCleanSamcommandforsamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createSortcommandforbamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createBAIfileforbamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createFLAGSTATcommandforfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createMERGEcommandformousebamfiles.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createBAIfileforbamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createDEDUPcommandsforfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createTARGETrealignmentintervallist.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createINDELREALIGNMENTcommands.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createPICARDQCcommandforbamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createFASTQCcommandforbamfile.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createCOVERAGEDEPTHcommands.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createMPILEUPfilesforbam.R')

# set path to files
pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/cell.xlsx"


#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/DN19266cellCSCE10_1f.xlsx"


#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/DN19306all.xlsx"

#pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"

mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)

#pathtosamfiles <- "/camp/lab/swantonc/working/bakkerb/DN19306/"
#mousedata <- read.xlsx("/camp/lab/swantonc/working/bakkerb/DN19306/Debbie_WES_mice.xlsx", sheetIndex = 1, stringsAsFactors=FALSE)




mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
paddedexonbedfile <- "/camp/lab/swantonc/working/bakkerb/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"
#paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/S0276129_Paddedliftovermm10.bed"


## specify NGSIDs
NGSIDs <- mousedata$Genomics.ID

for(i in 1:length(NGSIDs)){
  
  NGSID <- NGSIDs[i]
  mouse.ID  <- mousedata$Mouse.Name[mousedata$Genomics.ID%in%NGSID]
  sample.ID <- mousedata$Sample.Name[mousedata$Genomics.ID%in%NGSID]
  
  if(length(NGSID)<1 | length(mouse.ID)<1 | length(sample.ID)<1){
    print(paste("Mouse details not found for NGSID:", NGSID))
    next()
  }else{
    
    ###################################        REMOVE ALL UNNECESSARY INTERMEDIATE BAMs        ######################################### 
    
    mousepath <- paste0(pathtosamfiles, mouse.ID)
    intermediate.files <- list.files(mousepath, recursive = FALSE, full.names = T, pattern = NGSID)
    intermediate.files.to.delete <- intermediate.files[which(!grepl("processed", intermediate.files))]
    intermediate.files.to.delete <- intermediate.files.to.delete[!grepl(".sh", intermediate.files.to.delete)]
    
    intermediate.files.to.keep <- intermediate.files[which(grepl("processed", intermediate.files))]
    intermediate.files.to.keep <- unique(intermediate.files.to.keep)
    if(length(intermediate.files.to.keep)!=3){
      print(paste("Likely missing BAM/mpileup for NGSID/sampleID or already cleaned up:", NGSID, sample.ID))
      next()
    }
    if(length(intermediate.files.to.delete)>=1){
      rm.command <- paste("rm", intermediate.files.to.delete)
      print(rm.command)
      for(j in 1: length(rm.command)){
        system(rm.command[j])
      }
    }
    
    ###################################        RENAME THE BAM FILES TO SAMPLE IDS        #########################################
    new.intermediate.files.to.keep.names <- gsub(NGSID, sample.ID, intermediate.files.to.keep)
    mv.command <- paste("mv", intermediate.files.to.keep, new.intermediate.files.to.keep.names)
    print(mv.command)
    for(j in 1:length(mv.command)){
        system(mv.command[j])
    }

  }
}
  
  
  
  
