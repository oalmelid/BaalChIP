#BaalChIP: functions to compute alleleCounts
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


get_snp_ranges <- function(snps) {
    #suppressPackageStartupMessages(require(GenomicRanges))
    GRanges(snps$CHROM, IRanges(snps$POS, snps$POS,names=snps$ID), REF=snps$REF, ALT=snps$ALT)
}

bed2ranges <- function(bed) {
    #suppressPackageStartupMessages(require(GenomicRanges))
    gi.ranges <- GRanges(bed[,1], IRanges(as.numeric(bed[,2])+1, as.numeric(bed[,3])))
    gi.ranges
}

get_snps_in_GI <- function(snpfile, bedfile){
    #suppressPackageStartupMessages(require(GenomicRanges))

    snps <- read.delim(snpfile, stringsAsFactors=FALSE, header=TRUE)
    bed <- read.delim(bedfile, stringsAsFactors=FALSE, header=FALSE)
    snp.ranges <- get_snp_ranges(snps)
    gi.ranges <- bed2ranges(bed)
    ov <- suppressWarnings(overlapsAny(snp.ranges, gi.ranges, ignore.strand = TRUE))

    #Snps in genomic intervals
    sigi.ranges <- snp.ranges[ov,]
    sigi.ranges
}

filter_sigi <- function(snpfile, bedfile) {
    #get snps in bed file
    if (!is.null(bedfile)) {
        sigi.ranges <- get_snps_in_GI(snpfile, bedfile)
    }

    #if no genomic regions are given, get snp range object
    else {
        snps <- read.delim(snpfile, stringsAsFactors=FALSE, header=TRUE)
        sigi.ranges <- get_snp_ranges(snps)
    }

    sigi.ranges
}

tablenucs <- function(pileupres) {
    nucleotides <- levels(pileupres$nucleotide)
    res <- split(pileupres, pileupres$seqnames)
    res <- lapply(res, function (x) {split(x, x$pos)})
    res <- lapply(res, function (positionsplit) {
        nuctab <- lapply(positionsplit, function(each) {
                        chr = as.character(unique(each$seqnames))
                        pos = as.character(unique(each$pos))
                        tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
                        c(chr,pos, tablecounts)
                    })
        nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=FALSE)
        rownames(nuctab) <- NULL
        nuctab
    })
    res <- data.frame(do.call("rbind", res),stringsAsFactors=FALSE)
    rownames(res) <- NULL
    colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
    res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
    res
}


get_allele_counts <- function(bamfile, snp.ranges, returnRanges=FALSE,min_base_quality=10,min_mapq=15) {

    #suppressPackageStartupMessages(require(Rsamtools))
    #suppressPackageStartupMessages(require(GenomicRanges))

    #match sequences between snp.ranges and bamfile
    #snp.ranges cannot contain ranges that are not in the bamheader, otherwise pileup will give an error
    bf <- BamFile(bamfile) #create a bamfile instance
    bam_seqlengths <- seqlengths(bf)
    olap <- suppressWarnings(overlapsAny(snp.ranges, GRanges(names(bam_seqlengths),
                                         IRanges(rep(1, length(bam_seqlengths)), as.numeric(bam_seqlengths)))))

    if (sum(olap) == 0) {stop(c("SNPs do not overlap sequence names in bam file\nBAM seqlengths:\n",
    paste(names(bam_seqlengths),collapse=","),"\n\nSNP ranges:\n",paste(levels(seqnames(snp.ranges)),collapse=",")))}

    snp.ranges <- snp.ranges[olap]

    #compute pileup
    bf <- BamFile(bamfile) #create a bamfile instance
    param <- ScanBamParam(which=snp.ranges)
    pileupres <- pileup(bf, scanBamParam=param, pileupParam=PileupParam(
    max_depth=1000,
    min_mapq=min_mapq, min_base_quality=min_base_quality,
    distinguish_nucleotides=TRUE, distinguish_strands=FALSE, ignore_query_Ns=FALSE))

    #get table of nucleotide counts
    if (nrow(pileupres)==0) {return(NULL)} #there are no reads overlapping any snp.ranges
    nuctab <- tablenucs(pileupres)


    #merge with snp.ranges to get info about REF, ALT alleles
    snps <- data.frame("names"= names(snp.ranges), "seqnames"=seqnames(snp.ranges),
                               "start"=start(snp.ranges),
                               "REF"=as.character(values(snp.ranges)$REF),
                               "ALT"=as.character(values(snp.ranges)$ALT),
                               stringsAsFactors=FALSE)

    #snps <- as.data.frame(snp.ranges)
    #snps$names <- rownames(snps)
    snps <- merge(snps, nuctab, by=c("seqnames","start"))
    #snps$REF <- as.character(snps$REF)
    #snps$ALT <- as.character(snps$ALT)


    ref.counts <- sapply(seq_len(nrow(snps)), function(x) {snps[x,snps[x,"REF"]]})
    alt.counts <- sapply(seq_len(nrow(snps)), function(x) {snps[x,snps[x,"ALT"]]})
    total.counts <- ref.counts + alt.counts
    total.counts.withForeignReads <- rowSums(snps[, (colnames(nuctab)[-c(1,2)])])
    foreign.counts <-  total.counts.withForeignReads - total.counts

    #-------------------- return d.frame --------------------#
    allelecounts <- data.frame("ID"=snps$names, "CHROM"=snps$seqnames,
                               "POS"=snps$start,
                               "REF"=snps$REF, "ALT"=snps$ALT,
                               "REF.counts"=ref.counts,
                               "ALT.counts"=alt.counts,
                               "Total.counts"= total.counts,
                               "Foreign.counts"=foreign.counts,
                               "AR" = ref.counts / total.counts,
                               stringsAsFactors=FALSE)

    allelecounts <- allelecounts[allelecounts$Total.counts > 0 ,, drop=FALSE] #eliminate if it is all zero!
    if (!returnRanges) {return(allelecounts)}


    #-------------------- return ranges --------------------#
    if (returnRanges) {
        sigi.ranges <- makeGRangesFromDataFrame(allelecounts, seqnames.field="CHROM",
                                            start.field="POS",end.field="POS",
                                            keep.extra.columns = TRUE)
        names(sigi.ranges) <- allelecounts$ID
        return(sigi.ranges)
    }
}

applyAlleleCountsPerBam <- function(samples, hets, min_base_quality=min_base_quality, min_mapq=min_mapq, verbose=TRUE) {
    cells <- unique(samples[["group_name"]])
    res_per_bam <- list()
    res_per_bam <- lapply(cells, function(x) {res_per_bam[[x]] <- list()})
    names(res_per_bam) <- cells
    if (verbose) {
        message("-computing allele counts per BAM")
        pb <- txtProgressBar(min = 0, max = nrow(samples), style = 3)
    }
    for (rownr in seq_len(nrow(samples))) {

        x <- samples[rownr,,drop=FALSE]

        #get SNPs in genomic intervals (peaks, genes)
        snpfile = hets[[x[["group_name"]]]]
        sigi.ranges <- filter_sigi(snpfile=snpfile, bedfile = x[["bed_name"]])

        #Count frequency of Ref and alternative alleles
        sigi.ranges <- get_allele_counts(bamfile = x[["bam_name"]], snp.ranges = sigi.ranges, returnRanges=TRUE, min_base_quality=min_base_quality, min_mapq=min_mapq)
        res_per_bam[[x[["group_name"]]]][[x[["sampleID"]]]] <- list("sigi"=sigi.ranges)

        #set progress bar
        if (verbose) {setTxtProgressBar(pb, rownr)}
    }
    if (verbose) {close(pb)}

    return(res_per_bam)

}
