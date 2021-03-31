##SCRIPT_PLOT_MUT_SIG and heatmaps
# fail CSCE10_1c CSCE9_1e  CSCE4_1b CSCE9_1j CSCE10_2a
# heatmap production fails if there are fewer than 3 mutations
library(xlsx)
library(Heatplus)
library(gplots)
library(RColorBrewer)
library(ape)
library(phangorn)
library(deconstructSigs)
library('BSgenome.Mmusculus.UCSC.mm10')

source("/camp/project/proj-tracerx-lung/tracerx/_PIPELINE/tracerx-exome-pipeline/tracerx.functions.camp.R")

##Debbie
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
#mousedata <- read.xlsx("/camp/lab/swantonc/working/bakkerb/DN19306/Debbie_WES_mice.xlsx", sheetIndex = 1, stringsAsFactors=FALSE)

# references
paddedexonbedfile <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/twist_mouse_exome_targets_v1p1_mm10.bed"
mouserefpath <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/mm10/mm10.fa"

dbsnp <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/S0276129/mgp.v4.snps.dbSNP_chr.vcf"
blacklist.path <- "/camp/lab/swantonc/working/albakim/MouseImmunoediting/AnnotationFiles/mm10.blacklist.bed"
tmp.dir <- "/camp/lab/swantonc/working/bakkerb/DN19306/tmp/"
annovar.path <- '/camp/lab/swantonc/working/albakim/MouseImmunoediting/Annovar/table_annovar.pl'
annovar.params <- "-remove -buildver mm10 -protocol refGene,cytoBand,genomicSuperDups,bed,snp142 -operation g,r,r,r,f -bedfile mm10_blacklisted.bed /camp/lab/swantonc/working/albakim/MouseImmunoediting/Annovar/mm10db/ -nastring NA"

## specify mice
mice <- as.character(unique(mousedata$Mouse.Name))

for(i in 1:length(mice)){
        m <- mice[i]
#        print(m)
        
        ##check if SNV file exists - if not print error and move on 
        SNVdir <- paste0(pathtosamfiles, m, "/SNV")
        if(!dir.exists(SNVdir)){
                print(paste(m , "has no SNV folder"))
                next()
        }
        
        ##check if the SNV .xls output exists
        xlsfiles <- list.files(SNVdir, pattern = "SNV.xls", full.names = TRUE) # added because of custom xls files
        if(length(xlsfiles)<1){
                print(paste(m , "has missing .xls files"))
                next()
        }
        
        ##plot heatmaps - filtered and unfiltered
        snvFile <- xlsfiles[grepl("Exome.SNV.xls", xlsfiles)]
#         snvFile <- xlsfiles[grepl(".filtered.SNV.xls", xlsfiles)]
        if(any(!file.exists(snvFile))){
                print(paste0("missing varFiles for mouse: ", m))
                next()
	}

         cat(i, " ", snvFile, " ")

                if(grepl("filtered", snvFile)){
                        label <- "filtered"
                }else{
                        label <- "unfiltered"
                }
                mutations <- read.table(snvFile, sep="\t", header = T, stringsAsFactors = F)
                if (dim(mutations)[1] ==0 ) {
                  cat(" skip 0  \n")
                  next()
                }

                mutations <- mutations[mutations$is_SNV,]
                if (dim(mutations)[1] < 50 ) {
#                if (dim(mutations)[1] < 15 ) {
#Warning message:
#In mut.to.sigs.input(mut.ref = to.plot, sample.id = "sample", chr = "chr",  :
#  Some samples have fewer than 50 mutations:

                  cat(" skip < 50  \n")
                  next()
                }
                cat("plot ",snvFile, " ", dim(mutations)[1], "\n")
                mutations$sample.id <- "1"
                to.plot <- mutations[,c("sample.id", "chr", "start", "ref", "var")]
                colnames(to.plot) <- c("sample", "chr", "pos", "ref", "alt")
                sigs.input <- mut.to.sigs.input(mut.ref = to.plot,
                                                sample.id = "sample", 
                                                chr = "chr", 
                                                pos = "pos", 
                                                ref = "ref", 
                                                alt = "alt", 
                                                bsg=BSgenome.Mmusculus.UCSC.mm10)
                
                tumour<- NULL
                tumour$context <- colnames(sigs.input)
                tumour<- as.data.frame(tumour)
                tumour$fraction <- as.numeric(sigs.input[1,]/sum(sigs.input[1,]))
                tumour$trinuc <- substr(tumour$context,3,5)
                tumour$colour <- ifelse(tumour$trinuc=="C>A","#999999", NA )
                tumour$colour <- ifelse(tumour$trinuc=="C>G","#E69F00", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="C>T","#56B4E9", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>A","#009E73", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>C","#F0E442", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>G","#0072B2", tumour$colour  )
                FFPEtumour <- tumour
                pdf(file=paste0(SNVdir, "/", m, ".", label, ".mutSig.barplot.pdf"), height =4 , width=7);
                barplot(tumour$fraction, col=tumour$colour, names.arg=tumour$context, las=2, cex.names = 0.5, space = 0.5)
                dev.off()
     }





for(i in 1:length(mice)){
        m <- mice[i]
#        print(m)
        
        ##check if SNV file exists - if not print error and move on 
        SNVdir <- paste0(pathtosamfiles, m, "/SNV")
        if(!dir.exists(SNVdir)){
                print(paste(m , "has no SNV folder"))
                next()
        }
        
        ##check if the SNV .xls output exists
        xlsfiles <- list.files(SNVdir, pattern = "SNV.xls", full.names = TRUE) # added because of custom xls files
        if(length(xlsfiles)<1){
                print(paste(m , "has missing .xls files"))
                next()
        }
        
        ##plot heatmaps - filtered and unfiltered
#        snvFile <- xlsfiles[grepl("Exome.SNV.xls", xlsfiles)]
         snvFile <- xlsfiles[grepl(".filtered.SNV.xls", xlsfiles)]
        if(any(!file.exists(snvFile))){
                print(paste0("missing varFiles for mouse: ", m))
                next()
	}

         cat(i, " ", snvFile, " ")

                if(grepl("filtered", snvFile)){
                        label <- "filtered"
                }else{
                        label <- "unfiltered"
                }
                mutations <- read.table(snvFile, sep="\t", header = T, stringsAsFactors = F)
                if (dim(mutations)[1] ==0 ) {
                  cat(" skip 0  \n")
                  next()
                }

                mutations <- mutations[mutations$is_SNV,]
                if (dim(mutations)[1] < 50 ) {
#                if (dim(mutations)[1] < 15 ) {
#Warning message:
#In mut.to.sigs.input(mut.ref = to.plot, sample.id = "sample", chr = "chr",  :
#  Some samples have fewer than 50 mutations:

                  cat(" skip < 50  \n")
                  next()
                }
                cat("plot ",snvFile, " ", dim(mutations)[1], "\n")
                mutations$sample.id <- "1"
                to.plot <- mutations[,c("sample.id", "chr", "start", "ref", "var")]
                colnames(to.plot) <- c("sample", "chr", "pos", "ref", "alt")
                sigs.input <- mut.to.sigs.input(mut.ref = to.plot,
                                                sample.id = "sample", 
                                                chr = "chr", 
                                                pos = "pos", 
                                                ref = "ref", 
                                                alt = "alt", 
                                                bsg=BSgenome.Mmusculus.UCSC.mm10)
                
                tumour<- NULL
                tumour$context <- colnames(sigs.input)
                tumour<- as.data.frame(tumour)
                tumour$fraction <- as.numeric(sigs.input[1,]/sum(sigs.input[1,]))
                tumour$trinuc <- substr(tumour$context,3,5)
                tumour$colour <- ifelse(tumour$trinuc=="C>A","#999999", NA )
                tumour$colour <- ifelse(tumour$trinuc=="C>G","#E69F00", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="C>T","#56B4E9", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>A","#009E73", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>C","#F0E442", tumour$colour  )
                tumour$colour <- ifelse(tumour$trinuc=="T>G","#0072B2", tumour$colour  )
                FFPEtumour <- tumour
                pdf(file=paste0(SNVdir, "/", m, ".", label, ".mutSig.barplot.pdf"), height =4 , width=7);
                barplot(tumour$fraction, col=tumour$colour, names.arg=tumour$context, las=2, cex.names = 0.5, space = 0.5)
                dev.off()
     }

