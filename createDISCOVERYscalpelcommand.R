createDISCOVERYscalpelcommand <- function(NGSID, GERMLINE, tumour, germline, reference, bedfile, SCALPELdirectory){
        outputdirectory <- SCALPELdirectory
        command <- paste("/camp/lab/swantonc/working/naganoa/Application/myscalpel/scalpel-0.5.4/scalpel-discovery --somatic --normal ", germline, " --tumor ", tumour, " --bed ", bedfile, " --ref ", reference, " --dir ", outputdirectory, " --numprocs 4 ", sep="")
        command <- gsub(pattern = "/Volumes/lab-swantonc/", replacement = "/camp/lab/swantonc/", command)
        return(command)
}