#BaalChIP: filtering functions for QC
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

get_readlen <- function(bamfile, snp.ranges) {
    #suppressPackageStartupMessages(require(GenomicAlignments))
    #suppressPackageStartupMessages(require(Rsamtools))
    if (length(snp.ranges) >= 100) {Param <- ScanBamParam(which= snp.ranges[1:100])}
    if (length(snp.ranges) < 100)  {Param <- ScanBamParam(which= snp.ranges)}
    bf <- BamFile(bamfile) #create a bamfile instance
    temp <- readGAlignments(bf,param=Param)
    readlength=round(mean(qwidth(temp)))
    readlength
}

applyReadlenPerBam <- function(samples, res_per_bam, verbose=TRUE) {

    if (verbose) {
        message("-getting read lengths per sample")
        pb <- txtProgressBar(min = 0, max = nrow(samples), style = 3)
    }
    readlens <- lapply(seq_len(nrow(samples)), function(rownr) {
        #for (rownr in seq_len(nrow(samples))) {
        x <- samples[rownr,,drop=FALSE]

        x1 <- res_per_bam[[x[["group_name"]]]][[x[["sampleID"]]]] #get lastset
        lastset <- x1[[length(x1)]]
        sigi.ranges <- lastset

        #readlens <- c(readlens, get_readlen(bamfile=x[["bam_name"]], snp.ranges=sigi.ranges))
        result <- get_readlen(bamfile=x[["bam_name"]], snp.ranges=sigi.ranges)
        if (verbose) {setTxtProgressBar(pb, rownr)}
        return(result)
    })
    if (verbose) {close(pb)}
    readlens <- unlist(readlens)
    return(readlens)
}

get_lastset <- function(res_per_bam) {

    lastset <- lapply(names(res_per_bam) , function(group_name) {
        #for (group_name in names(res_per_bam)) {
        counts_per_group <- res_per_bam[[group_name]]
        l <- lapply(counts_per_group, function(x) {x[[length(x)]]})
        l <- lapply(l , function(x) {
                    data.frame("ID"=names(x),
                               "seqnames"=seqnames(x),
                               "start"=start(x),
                               "REF"=values(x)$REF,
                               "ALT"=values(x)$ALT)
        })
        l <- lapply(l, subset, select=c("ID","seqnames","start","REF","ALT"), drop=FALSE)
        l <- lapply(l, function(x) {rownames(x) <- NULL; x})
        #lastset[[group_name]] <- l
        return(l)
    })
    names(lastset) <- names(res_per_bam)

    lastset <- lapply(lastset, function(l) { l <- do.call("rbind",l); rownames(l) <- NULL; l})
    lastset <- unique(do.call("rbind",lastset))
    rownames(lastset) <- NULL
    colnames(lastset) <- c("ID","CHROM","POS","REF","ALT")
    return(lastset)
}

makeARGSTRING <- function(alignmentSimulArgs) {
    if (is.null(alignmentSimulArgs)) {return("")}
    names(alignmentSimulArgs) <- NULL
    paste(alignmentSimulArgs, collapse=" ")
}

run_sim <- function (snpframe, readlenvector, outputPath, simulation_script, alignmentSimulArgs=NULL, verbose=TRUE) {

    #output snp file name
    snpout <- paste0(outputPath,"_snplist.txt")

    #save snp.ranges into a file to be used by perl scripts
    a <- cbind(snpframe[,c("ID","CHROM","POS"),drop=FALSE], "+",snpframe[,c("REF","ALT"),drop=FALSE])
    write.table(a, file = snpout, sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

    #run external simulation script
    if (!is.null(simulation_script)) {

        if (verbose) {
            message("-running simulations for each read length")
            pb <- txtProgressBar(min = 0, max = length(readlenvector), style = 3)
        }

        for (i in seq_len(length(readlenvector))) {


            #imput arguments of run_simulations.sh
            readlen <- readlenvector[[i]]
            fastaout <- paste0(outputPath,"_aln",readlen,".fasta")
            bamout <- paste0(outputPath,"_aln",readlen,".bam")

            perlScript <- file.path(system.file("extra", package="BaalChIP"),"get.overlaps.v2_chrXY.perl")

            ARGSTRING <- makeARGSTRING(alignmentSimulArgs)

            command_to_run <- paste(simulation_script,perlScript,readlen,snpout,fastaout,bamout,ARGSTRING)

            r <- system(command_to_run, intern=TRUE)
            if (!is.null(attr(r,"status"))) {
                stop("error running external script")
            }
            #debug: if (r != 0) {stop("error running simulation script")}
            #debug: "bsub -q fmlab -R rusage[mem=5000]"

            if (verbose) {setTxtProgressBar(pb, i)}

        }
        if (verbose) {close(pb)}
    }
}

get_simcounts <- function(snpframe, readlenvector, outputPath, verbose=TRUE) {

    #get simcounts
    if(verbose) {
        message("-calculating allelecounts for simulated reads")
        pb <- txtProgressBar(min = 0, max = length(readlenvector), style = 3)
    }

    simcounts <-  lapply(seq_len(length(readlenvector)), function(i) {
        r <- readlenvector[[i]]
        bamout <- paste0(outputPath,"_aln",r,".bam")
        res <- get_allele_counts(bamfile = bamout, snp.ranges = get_snp_ranges(snpframe), returnRanges=FALSE, min_base_quality=0, min_mapq=0)
        if (verbose) {setTxtProgressBar(pb, i)}
        return(res)
    })
    if (verbose) {close(pb)}

    names(simcounts) <- readlenvector
    simcounts
}


applySim <- function(samples, res_per_bam, simul_output, simulation_script= "local", alignmentSimulArgs=NULL, testingScript=FALSE, verbose=TRUE) {

     if (simulation_script == "local") {simulation_script = system.file("extra/run_simulations.sh",package="BaalChIP")}
     lastset <- get_lastset(res_per_bam) #get the data frame for all unique SNPs from the last round of filters
     readlenvector <- unique(samples$readlen) #vector of read lens

     message("-detected readlengths", readlenvector, "\n")

     #run_sum
     if (!testingScript) { #will save fasta+bam files in simul_output
        run_sim(snpframe=lastset, readlenvector, outputPath=simul_output, simulation_script=simulation_script, alignmentSimulArgs=alignmentSimulArgs, verbose=verbose)
     }

     simcounts <- get_simcounts(snpframe=lastset, readlenvector, outputPath=simul_output, verbose=verbose)

     #simulation_stats
     total_nr_snps <- nrow(lastset)
     simulation_stats <- lapply(names(simcounts), function(r) {
                                                sim <- simcounts[[r]]
                                                r <- as.numeric(r)
                                                perc_right = sum(sim$AR==0.5 & sim$REF.counts == r*2 & sim$ALT.counts == r*2)
                                                c("readslen"=r,"perc_right"=100*perc_right / total_nr_snps)
                                                })


      simulation_stats <- data.frame(do.call("rbind", simulation_stats), stringsAsFactors=FALSE)

      return(list("simcounts" = simcounts, "simulation_stats" = simulation_stats))
}

filter_intbias <- function(snp.ranges, sim, r){

     #only keep SNPs that have the right number of reads mapped to them
     keep <- as.character(sim[sim$AR==0.5 & sim$REF.counts == r*2 & sim$ALT.counts == r*2, "ID"])

     #filter out SNPs that are not 50-50 in the simulated dataset
     snp.ranges[names(snp.ranges) %in% keep,]

}

applyIntBiasFilterPerBam <- function(samples, res_per_bam, simcounts, verbose=TRUE) {
    if (verbose) {
        message("-filtering intrinsic bias")
        pb <- txtProgressBar(min = 0, max = nrow(samples), style = 3)
    }
    for (rownr in seq_len(nrow(samples))) {
            x <- samples[rownr,,drop=FALSE]

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
            if (verbose) {setTxtProgressBar(pb, rownr)}
        }
        if (verbose) {close(pb)}
        return(res_per_bam)
}



