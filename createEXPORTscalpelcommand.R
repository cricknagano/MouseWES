createEXPORTscalpelcommand <- function(NGSID, GERMLINE, tumour, germline, reference, bedfile, SCALPELdirectory){
        db <- paste(SCALPELdirectory, "/main/somatic.db", sep="")
        outputfile <- paste(SCALPELdirectory, "/", NGSID,"vs", GERMLINE, "_scalpelexport.vcf", sep="" )
        command <- paste("/camp/lab/swantonc/working/naganoa/Application/myscalpel/scalpel-0.5.4/scalpel-export --somatic --db ", db, " --bed ", bedfile, " --ref ", reference, " --min-vaf-tumor 0.001 > ", outputfile, sep="")
        command <- gsub(pattern = "/Volumes/lab-swantonc/", replacement = "/camp/lab/swantonc/", command)
        return(command)
}