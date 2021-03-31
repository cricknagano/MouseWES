## SOFT LINK THE CNV OUTPUT TO CORRECT FOLDERS
library(xlsx)

source('/camp/lab/swantonc/working/albakim/MousePipeline/createCNVkitcommand.R')
source('/camp/lab/swantonc/working/albakim/MousePipeline/createCNVkitcommandrefcnn.R')

## SET ALL PATHS
pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/"
mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/cell.xlsx"

#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19266/DN19266cell.xlsx"

#pathtosamfiles<- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_3/DN19306all.xlsx"
#pathtosamfiles <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/"
#mousedatapath <- "/camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN19306_2/Debbie_WES_info.xlsx"

mousedata <- read.xlsx(mousedatapath, sheetIndex = 1, stringsAsFactors=FALSE)

#pathtosamfiles <- "/camp/lab/swantonc/working/bakkerb/DN19306/"
#mousedata <- read.xlsx("/camp/lab/swantonc/working/bakkerb/DN19306/Debbie_WES_mice.xlsx", sheetIndex = 1)


mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"
dbsnp <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/mgp.v4.snps.dbSNP_chr.vcf"
path <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/FASTQfolder/"
mybaits <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/cnvkit/my_UCSCtargets.bed"
accessfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/cnvkit/access-10kb.mm10.bed"
annotatefile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10_refFlat.txt"

NGSIDs <- as.character(mousedata$Genomics.ID)
mice <- as.character(unique(mousedata$Mouse.Name))
tumour.sample.names <- as.character(unique(mousedata$Sample.Name[mousedata$sample.type%in%"tumour"]))

CNVfolder <- paste0(pathtosamfiles, "CNV")

for(i in 1:length(mice)){
  m <- mice[i]
  f <- list.files(CNVfolder, pattern = m, full.names = TRUE, recursive = F)
  new.path <- paste0(pathtosamfiles, m, "/CNV")
  if(!file.exists(new.path)){
    dir.create(new.path)
  }
  for(j in 1:length(f)){
    print(f[j])
    system(paste("ln -s", f[j], new.path))
  }
}
