##SCRIPT_FOR_CNV_ANALYSIS
library(xlsx)

source('/camp/lab/swantonc/working/albakim/MousePipeline/createCNVkitcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createCNVkitcommandrefcnn.R')

## SET ALL PATHS
pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/DN19306all.xlsx"


#pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"
mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)

#pathtosamfiles <- "/camp/lab/swantonc/working/bakkerb/DN19306/"
#mousedata <- read.xlsx("/camp/lab/swantonc/working/bakkerb/DN19306/Debbie_WES_mice.xlsx", sheetIndex = 1)
#mousedata <- mousedata[!grepl("TOP UP", mousedata$Comments),]


mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"
#paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/S0276129_Paddedliftovermm10.bed"

dbsnp <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/mgp.v4.snps.dbSNP_chr.vcf"
path <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/FASTQfolder/"
mybaits <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/cnvkit/my_UCSCtargets.bed"
accessfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/cnvkit/access-10kb.mm10.bed"
annotatefile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10_refFlat.txt"

NGSIDs <- as.character(mousedata$Genomics.ID)
mice <- as.character(unique(mousedata$Mouse.Name))
tumour.sample.names <- as.character(unique(mousedata$Sample.Name[mousedata$sample.type%in%"tumour"]))


#########################################################################
##create supermouse shell script (ARM LEVEL)
#########################################################################
allbams <- list.files(paste0(pathtosamfiles, unique(mousedata$Mouse.Name)), recursive = FALSE, full.names = T, pattern = ".bam$")
allGLdirs <- allbams[grepl("Tail", allbams)]
allTUMdirs <- allbams[!allbams%in%allGLdirs]

## Generate script
shellDir <- paste0(pathtosamfiles, "/batchednorm.sh")
numberofcpus <- 8
erroroutput <- paste0(pathtosamfiles, "/", Sys.Date(), ".CNVkit.batchednorm.error.txt")

cat("#!/bin/bash", file = shellDir, sep='\n\n', append = FALSE, fill = FALSE)

cat(paste("#SBATCH --time=3-00:00:0"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
cat(paste("#SBATCH -e ", erroroutput, sep=""), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
##cat(paste("#SBATCH --exclusive", sep=""), file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)
cat(paste("#SBATCH --partition=cpu"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)
cat(paste("#SBATCH --cpus-per-task=8"), file = shellDir, sep='\n', append = TRUE, fill = FALSE)

cat("module load cnvkit/0.9.5-foss-2016b-Python-2.7.12", file = shellDir, sep='\n\n', append = TRUE, fill = FALSE)



  outputdir <- paste0(pathtosamfiles, "CNV/")
  if(!dir.exists(outputdir)){
    dir.create(outputdir)
  }
  cat(paste("srun -c ", numberofcpus, " -n 1 -J cnvkit time cnvkit.py batch ", paste(allTUMdirs, collapse=" "), " --normal ", paste(allGLdirs, collapse = " "), 
            " --targets ", paddedexonbedfile,  " --annotate /camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10_refFlat.txt --fasta /camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa --access /camp/lab/swantonc/working/albakim/MouseImmunoediting/cnvkit/access-10kb.mm10.bed", 
            " --output-reference ", paste0(outputdir,"/batched.reference"),
            " --output-dir ", outputdir,
            " -p 8",
            sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
  

cnr.paths <- paste0(outputdir, tumour.sample.names, "_processed.cnr")
cns.paths <- paste0(outputdir, tumour.sample.names, "_processed.cns")
cns.ed.path <- paste0(outputdir, tumour.sample.names, "_processed.ed.cns")
scatter.path <- paste0(outputdir, tumour.sample.names, "_processed.scatter.pdf")

for( i in 1:length(tumour.sample.names)){
  cat(paste("srun -c 1 -n 1 -J cnvseg time cnvkit.py segment ", cnr.paths[i]," -o ", cns.ed.path[i] ," -m none --drop-low-coverage", sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
  cat(paste("srun -c 1 -n 1 -J cnvsca time cnvkit.py scatter ", cnr.paths[i], " -s ", cns.ed.path[i], " --y-min -3 --y-max 3 -o ", scatter.path[i], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
  cat(paste("srun -c 1 -n 1 -J cnvdia time cnvkit.py diagram -s ", cns.ed.path[i], " ", cnr.paths[i], sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)
}

all.cns.paths <- paste0(outputdir, "*.ed.cns")
heatmap.pdf.path <- paste0(outputdir, "concatenated.heatmap.pdf")
cat(paste("srun -c 1 -n 1 -J cnvhea time cnvkit.py heatmap ", all.cns.paths, " -o ", heatmap.pdf.path, sep=""), file=shellDir, sep='\n\n', append = TRUE, fill = FALSE)

cat("wait", file = shellDir, append = TRUE, fill = FALSE)


#########################################################################
##create supermouse shell script (kb-level)
#########################################################################
