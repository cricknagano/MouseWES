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
        print(m)
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
#        snvFile <- "output/DN19306_2/CSCE10_1c/SNV/CSCE10_1c.Exome.SNV.xls"
        if(any(!file.exists(snvFile))){
                print(paste0("missing varFiles for mouse: ", m))
                next()
	}
        dBase      <- read.delim(snvFile,header=T,sep="\t",stringsAsFactors=F);
        dBase      <- dBase[apply(dBase[,grep("VAF",colnames(dBase))],1,function(x){!any(is.na(x))}),];


        dBase.filt <- dBase[dBase$Use.For.Plots | dBase$Use.For.Plots.Indel,]
        if (dim(dBase.filt)[1] < 1) {
            cat("skip i=", i, " ", mice[i], "\n")
            next()
        }
 
        heatmapFilter <- sapply(strsplit(dBase.filt$Func.refGene,";"),function(x){ any(x %in% c("exonic","splicing"))}) &
                # exonic and/or splicing
                dBase.filt$ExonicFunc.refGene %in% 
         c("nonsynonymous SNV","stopgain SNV","stoploss SNV", "nonframeshift substitution", 
        "frameshift insertion", "frameshift substitution", "nonsynonymous", "unknown")



        


        clustDat.unfilt.VAF           <- dBase[,grep("VAF",colnames(dBase))];
        clustDat.unfilt.VAF           <- clustDat.unfilt.VAF[,!colnames(clustDat.unfilt.VAF) %in% "max.VAF"];
        # added 3-3-2020 (duplicate rownames not allowed)
        clustDat.unfilt.VAF.rn        <- paste(dBase$Gene.refGene," ",dBase$key,sep='');
        clustDat.unfilt.VAF           <- clustDat.unfilt.VAF[!duplicated(clustDat.unfilt.VAF.rn), ]
        clustDat.unfilt.VAF.rn        <- clustDat.unfilt.VAF.rn[!duplicated(clustDat.unfilt.VAF.rn)]
        rownames(clustDat.unfilt.VAF) <- clustDat.unfilt.VAF.rn
        colnames(clustDat.unfilt.VAF) <- sub(".VAF","",colnames(clustDat.unfilt.VAF));
        clustDat.unfilt.VAF.hcRow     <- hclust(dist(clustDat.unfilt.VAF,method="euclidean"),method="average");
        clustDat.unfilt.VAF.hcCol     <- hclust(dist(t(clustDat.unfilt.VAF),method="euclidean"),method="average");
        

        clustDat.unfilt.binary        <- clustDat.unfilt.VAF;
        clustDat.unfilt.binary[clustDat.unfilt.binary<1]  <- 0;
        clustDat.unfilt.binary[clustDat.unfilt.binary>=1] <- 1;
        #clustDat.unfilt.binary.hcRow  <- hclust(dist(clustDat.unfilt.binary,method="binary"),method="average");
        clustDat.unfilt.binary.hcCol  <- hclust(dist(t(clustDat.unfilt.binary),method="binary"),method="average");


#filter
        if ( dim(dBase.filt)[1] < 2 ) {
           cat("skip i=", i, " ", mice[i], "\n")
           next()
        }
        cat("plot i=", i, " ",  mice[i], " ", dim(dBase.filt)[1], " plotn=")

        dBase.filt.heat <-dBase.filt[heatmapFilter,]

        clustDat.filt.VAF             <- dBase.filt[,grep("VAF",colnames(dBase.filt))];
        clustDat.filt.VAF             <- clustDat.filt.VAF[,!colnames(clustDat.filt.VAF) %in% "max.VAF"];
        rownames(clustDat.filt.VAF)   <- paste(dBase.filt$Gene.refGene," ",dBase.filt$key,sep='');
        colnames(clustDat.filt.VAF)   <- sub(".VAF","",colnames(clustDat.filt.VAF));
        clustDat.filt.VAF.hcRow       <- hclust(dist(clustDat.filt.VAF,method="euclidean"),method="average");
        clustDat.filt.VAF.hcCol       <- hclust(dist(t(clustDat.filt.VAF),method="euclidean"),method="average");
        
        
        clustDat.filt.binary          <- clustDat.filt.VAF;
        clustDat.filt.binary[clustDat.filt.binary<1]  <- 0;
        clustDat.filt.binary[clustDat.filt.binary>=1] <- 1;
        #clustDat.filt.binary.hcRow    <- hclust(dist(clustDat.filt.binary,method="binary"),method="average");
        clustDat.filt.binary.hcCol    <- hclust(dist(t(clustDat.filt.binary),method="binary"),method="average");
        
        clustDat.heat.filt.VAF             <- dBase.filt.heat[,grep("VAF",colnames(dBase.filt))];
        clustDat.heat.filt.VAF             <- clustDat.heat.filt.VAF[,!colnames(clustDat.heat.filt.VAF) %in% "max.VAF"];
        rownames(clustDat.heat.filt.VAF)   <- paste(dBase.filt.heat$Gene.refGene," ",dBase.filt.heat$key,sep='');
        colnames(clustDat.heat.filt.VAF)   <- sub(".VAF","",colnames(clustDat.heat.filt.VAF));
        clustDat.heat.filt.VAF.hcRow       <- hclust(dist(clustDat.heat.filt.VAF,method="euclidean"),method="average");
        clustDat.heat.filt.VAF.hcCol       <- hclust(dist(t(clustDat.heat.filt.VAF),method="euclidean"),method="average");
        
        clustDat.heat.filt.binary          <- clustDat.heat.filt.VAF;
        clustDat.heat.filt.binary[clustDat.heat.filt.binary<1]  <- 0;
        clustDat.heat.filt.binary[clustDat.heat.filt.binary>=1] <- 1;
        clustDat.heat.filt.binary.hcCol    <- hclust(dist(t(clustDat.heat.filt.binary),method="binary"),method="average");
        
        
        # reorder the binary matrices
        
        clustDat.unfilt.binary <- reorderBinaryMatrix(clustDat.unfilt.binary);
        clustDat.filt.binary   <- reorderBinaryMatrix(clustDat.filt.binary);
        
        clustDat.heat.filt.binary   <- reorderBinaryMatrix(clustDat.heat.filt.binary);
        
        upgmaTree<- function(data=data,filename=filename) {
                phyDat.dat<-phyDat(data,type="USER",levels=c(0,1))
                treeDat<-upgma(dist(t(data),method="binary"));
                treeDat<-acctran(treeDat,phyDat.dat)
                #pdf(file=paste(filename,".pdf",sep=""))
                #plot.phylo(treeDat,no.margin=TRUE,lab4ut="horizontal",edge.width=3,cex=0.7,use.edge.length=TRUE,type="unrooted")
                #dev.off()
                #png(file=paste(filename,".png",sep=""))
                #plot.phylo(treeDat,no.margin=TRUE,lab4ut="horizontal",edge.width=3,cex=0.7,use.edge.length=TRUE,type="unrooted")
                #dev.off()
        }
        upgmaTree(data=as.data.frame(clustDat.unfilt.binary),filename=clustDat.unfilt.binary.unrooted.file)
        upgmaTree(data=as.data.frame(clustDat.filt.binary),filename=clustDat.filt.binary.unrooted.file)
        
        
        lineHeight         <- c(1,round(nrow(clustDat.heat.filt.VAF)/13));
        heatMargins        <- c(20,20);
        heatPlotHeight.png <- nrow(clustDat.heat.filt.VAF)*20;
        heatPlotWidth.png  <- 200*ncol(clustDat.heat.filt.VAF);
        heatPlotHeight.png[heatPlotHeight.png < 1000] <- 1000;
        heatPlotWidth.png[heatPlotWidth.png < 800]    <- 600;
        heatPlotHeight.pdf <- nrow(clustDat.heat.filt.VAF)*20/480*7;
        heatPlotWidth.pdf  <- 200*ncol(clustDat.heat.filt.VAF)/480*7;
        
        if (heatPlotHeight.png == 1000 | heatPlotWidth.png == 600) {
                lineHeight <- c(1,round(nrow(clustDat.heat.filt.VAF)/6));
                heatPlotHeight.pdf <- (heatPlotHeight.png / 480) * 7;
                heatPlotWidth.pdf  <- (heatPlotWidth.png / 480) * 7;
                heatMargins        <- c(20,20);
        }
        
        # include a conditional for the heatmap - it fails if there are fewer than 3 mutations
        cat(nrow(clustDat.heat.filt.binary), "\n")
        if(nrow(clustDat.heat.filt.binary) > 3) {
                
                rowColours.binary.filt <- rep("white",nrow(clustDat.heat.filt.binary));
                rowColours.binary.filt[dBase.filt.heat$driverCategory[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == "1A"] <- "black";
                rowColours.binary.filt[dBase.filt.heat$driverCategory[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == 1] <- "navy";
                rowColours.binary.filt[dBase.filt.heat$driverCategory[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == "2A"] <- "blue";
                rowColours.binary.filt[dBase.filt.heat$driverCategory[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == 2] <- "dodgerblue";
                rowColours.binary.filt[dBase.filt.heat$driverCategory[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == 3] <- "skyblue";
                rowColours.binary.filt[dBase.filt.heat$indelLungDriver[match(sub(".* ","",rownames(clustDat.heat.filt.binary)),dBase.filt.heat$key)] == "TRUE"] <- "red"; 
#                png(file=paste0(SNVdir, "/", m, ".filtered.binary.heatmap.png"),height=heatPlotHeight.png,width=heatPlotWidth.png);
#                heatmap.2(
#                        x               = as.matrix(clustDat.heat.filt.binary),
#                        scale           = "none",
#                        Rowv            = FALSE,
#                        Colv            = as.dendrogram(clustDat.heat.filt.binary.hcCol),
#                        key             = FALSE,
#                        col             = c("grey","darkred"),
#                        na.color        = "darkgrey",
#                        dendrogram      = "column",
#                        trace           = "none",
#                        tracecol        = "black",
#                        linecol         = NA,
#                        margins         = heatMargins,
#                        density.info    = "density",
#                        cexRow          = 1.5,
#                        cexCol          = 1.5,
#                        colsep          = 1:ncol(clustDat.heat.filt.binary),
#                        rowsep          = 1:nrow(clustDat.heat.filt.binary),
#                        RowSideColors   = rowColours.binary.filt,
#                        sepcolor        = "white",
#                        sepwidth        = c(0.01,0.01),
#                        lwid            = c(7,10),
#                        lhei            = lineHeight);
#                dev.off();
                
                pdf(file=paste0(SNVdir, "/", m, ".filtered.binary.heatmap.pdf"),height=heatPlotHeight.pdf,width=heatPlotWidth.pdf);
                heatmap.2(
                        x               = as.matrix(clustDat.heat.filt.binary),
                        scale           = "none",
                        Rowv            = FALSE,
                        Colv            = as.dendrogram(clustDat.heat.filt.binary.hcCol),
                        key             = FALSE,
                        col             = c("grey","darkred"),
                        na.color        = "darkgrey",
                        dendrogram      = "column",
                        trace           = "none",
                        tracecol        = "black",
                        linecol         = NA,
                        margins         = heatMargins,
                        density.info    = "density",
                        cexRow          = 1.5,
                        cexCol          = 1.5,
                        colsep          = 1:ncol(clustDat.heat.filt.binary),
                        rowsep          = 1:nrow(clustDat.heat.filt.binary),
                        RowSideColors   = rowColours.binary.filt,
                        sepcolor        = "white",
                        sepwidth        = c(0.01,0.01),
                        lwid            = c(7,10),
                        lhei            = lineHeight);
                dev.off();
                
        }
        
}





