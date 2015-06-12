#BaalChIP: filtering functions for QC
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
    
    
get_readlen <- function(bamfile, snp.ranges) {
    suppressPackageStartupMessages(require(GenomicAlignments))
    suppressPackageStartupMessages(require(Rsamtools))
    if (length(snp.ranges) >= 100) {Param <- ScanBamParam(which= snp.ranges[1:100])}
    if (length(snp.ranges) < 100)  {Param <- ScanBamParam(which= snp.ranges)}
    bf <- BamFile(bamfile) #create a bamfile instance 
    temp <- readGAlignmentsFromBam(bf,param=Param)
    readlength=round(mean(qwidth(temp)))
    readlength
}

applyReadlenPerBam <- function(samples, res_per_bam) {
    readlens <- c()
    cat("getting read lengths per sample\n\n")
    pb <- txtProgressBar(min = 0, max = nrow(samples), style = 3)
    for (rownr in 1:nrow(samples)) {
        x <- samples[rownr,]
        sigi.ranges <- res_per_bam[[x[["group_name"]]]][[x[["sampleID"]]]][["sigi"]]
        readlens <- c(readlens, get_readlen(bamfile=x[["bam_name"]], snp.ranges=sigi.ranges))
        setTxtProgressBar(pb, rownr)
    }
    close(pb)
    return(readlens)
}

get_lastset <- function(res_per_bam) {
    lastset <- lapply(res_per_bam, lapply, function (x) { x <- as.data.frame(x[[length(x)]])[,c("ID","seqnames","start","REF","ALT")]; rownames(x) <- NULL; x })
    lastset <- lapply(lastset, function(l) { l <- do.call("rbind",l); rownames(l) <- NULL; l})
    lastset <- unique(do.call("rbind",lastset))
    rownames(lastset) <- NULL
    colnames(lastset) <- c("ID","CHROM","POS","REF","ALT")
    lastset
}

run_sim <- function (snpframe, readlenvector, outputPath, simulation_script) {
    
    #output snp file name
    snpout <- paste0(outputPath,"_snplist.txt")
    
    #save snp.ranges into a file to be used by perl scripts
    a <- cbind(snpframe[,c("ID","CHROM","POS")], "+",snpframe[,c("REF","ALT")])
    write.table(a, file = snpout, sep=" ", row.names=F, col.names=F, quote=F)
    
    #run external simulation script
    if (!is.null(simulation_script)) {
        for (r in readlenvector) {
            fastaout <- paste0(outputPath,"_aln",r,".fasta")
            bamout <- paste0(outputPath,"_aln",r,".bam")
            command_to_run <- paste(simulation_script,r,snpout,fastaout,bamout)
            #print(command_to_run) #debug!
            r <- system(command_to_run, intern=FALSE)
            #if (r != 0) {stop("error running simulation script")}
            #"bsub -q fmlab -R rusage[mem=5000]" 
        }    
    }
}

get_simcounts <- function(snpframe, readlenvector, outputPath) {
    
    cat("calculating allelecounts for simulated reads\n\n")
    
    #get simcounts
    pb <- txtProgressBar(min = 0, max = length(readlenvector), style = 3)
    
    simcounts <-  lapply(1:length(readlenvector), function(i) {
        r <- readlenvector[[i]]
        bamout <- paste0(outputPath,"_aln",r,".bam")
        res <- get_allele_counts(bamfile = bamout, snp.ranges = get_snp_ranges(snpframe), returnRanges=FALSE, min_base_quality=0, min_mapq=0)    
        setTxtProgressBar(pb, i) #set progress bar
        return(res)
    })
    close(pb)
    
    names(simcounts) <- readlenvector
    simcounts
}

applySim <- function(samples, res_per_bam, simul_output, simulation_script= NULL, testingScript=FALSE) {
    
     cat ("running simulations...\n\n")
     lastset <- get_lastset(res_per_bam) #get the data frame for all unique SNPs from the last round of filters
     readlenvector <- unique(samples$readlen) #vector of read lens
     
     #run_sum
     if (!testingScript) {run_sim(snpframe=lastset, readlenvector, outputPath=simul_output, simulation_script=simulation_script)} #will save fasta+bam files in simul_output
     simcounts <- get_simcounts(snpframe=lastset, readlenvector, outputPath=simul_output) 
    
     #simulation_stats
     total_nr_snps <- nrow(lastset)
     simulation_stats <- lapply(names(simcounts), function(r) {
                                                                    sim <- simcounts[[r]]
                                                                    r <- as.numeric(r)
                                                                    perc_right = sum(sim$AR==0.5 & sim$REF.counts == r*2 & sim$ALT.counts == r*2)
                                                                    c("readslen"=r,"perc_right"=100*perc_right / total_nr_snps)
                                                                  })
        
        
      simulation_stats <- data.frame(do.call("rbind", simulation_stats))
      
      return(list("simcounts" = simcounts, "simulation_stats" = simulation_stats))
}

filter_intbias <- function(snp.ranges, sim, r){
     
     #only keep SNPs that have the right number of reads mapped to them
     keep <- as.character(sim[sim$AR==0.5 & sim$REF.counts == r*2 & sim$ALT.counts == r*2, "ID"]) 

     #filter out SNPs that are not 50-50 in the simulated dataset
     snp.ranges[names(snp.ranges) %in% keep,]
        
}

applyIntBiasFilterPerBam <- function(samples, res_per_bam, simcounts) {
    cat("filtering intrinsic bias\n\n")
    pb <- txtProgressBar(min = 0, max = nrow(samples), style = 3)
    for (rownr in 1:nrow(samples)) {
            x <- samples[rownr,]
            
            #get appropriate simulation data set (based on read length)
            sim <- simcounts[[as.character(x[["readlen"]])]]
            
            #get snps to filter out further
            sampledata <- res_per_bam[[x[["group_name"]]]][[x[["sampleID"]]]]
            snp.ranges <- sampledata[[length(sampledata)]]
            
            #filter intrinsic bias
            res <- filter_intbias(snp.ranges, sim, r=x[["readlen"]])
            
            #append the results
            res_per_bam[[x[["group_name"]]]][[x[["sampleID"]]]][["intbias"]] <- res
            
            #set progress bar
            setTxtProgressBar(pb, rownr)
        }      
        close(pb)
        return(res_per_bam)
} 



