#BaalChIP: filtering functions for QC
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

renameSeqlevels.allowNAs <- function(region.ranges, newStyle) {

    #filter only those coordinates that had a matching seqname
    newStyle <- newStyle[!is.na(newStyle)]
    region.ranges <- region.ranges[which(seqnames(region.ranges) %in% names(newStyle))]
    region.ranges <- GRanges(as.character(seqnames(region.ranges)), IRanges(start(region.ranges),end(region.ranges)))

    #map again
    newStyle <- mapSeqlevels(seqlevels(region.ranges), "NCBI")
    region.ranges <- renameSeqlevels(region.ranges, newStyle)
    region.seqlevels <- seqlevels(region.ranges)

    region.ranges
}

filter_regions <- function(snp.ranges, region.ranges, type=c("filterOut_overlapping_sites","keep_overlapping_sites")) {

    #check if seqlevels intersect
    snp.seqlevels <- seqlevels(snp.ranges)
    region.seqlevels <- seqlevels(region.ranges)

    #rename if no intersect
    if (length(intersect(snp.seqlevels, region.seqlevels))==0) {

        #rename region.ranges
        newStyle <- mapSeqlevels(region.seqlevels, "NCBI")
        region.ranges <- renameSeqlevels.allowNAs(region.ranges, newStyle)
        region.seqlevels <- seqlevels(region.ranges)

        #check again if intersect
        if (length(intersect(snp.seqlevels, region.seqlevels))==0) {
            warning("In intersect SNPs with Filtering regions:
            The 2 combined objects have no sequence levels in common.")
        }else{
            #message("-automatically matched between UCSC/NCBI style to find overlapping filtering regions")
        }
    }


    ov <- suppressWarnings(overlapsAny(snp.ranges,region.ranges,ignore.strand = TRUE))
    if (type == "filterOut_overlapping_sites") {return(snp.ranges[!ov,])}
    if (type == "keep_overlapping_sites") {return(snp.ranges[ov,])}
}

filter_genomicRegions <- function (snp.ranges, Regions, type=c("filterOut_overlapping_sites","keep_overlapping_sites")) {
    #loops through RegionsToFilter and returns a list of filtered SNPs
    #(applies hierarchy within RegionsToFilter)

    if (is.null(Regions)) {
        return(list(snp.ranges)) #nothing to filter

    }

    sigi.filtered <- list()
    snps2filter <- snp.ranges
    for (region_name in names(Regions)) {
            region <- Regions[[region_name]]
            filtered_snps <- filter_regions(snp.ranges = snps2filter, region.ranges = region, type=type)
            sigi.filtered[[region_name]] <- filtered_snps
            snps2filter = filtered_snps
    }
    return(sigi.filtered)
}


applyFiltersPerBam <- function(counts_per_bam, RegionsToFilter, RegionsToKeep, verbose=TRUE) {

    res_per_bam <- list()

    N <- sum(sapply(counts_per_bam, length))

    if (verbose) {
        message("-applying filters per BAM")
        pb <- txtProgressBar(min = 0, max = N, style = 3)
    }
    i=1

    for (group_name in names(counts_per_bam)) {

        for (sampleID in names(counts_per_bam[[group_name]])) {

            #GRanges of variants
            x <- counts_per_bam[[group_name]][[sampleID]] #get lastset
            sigi.ranges <- x[[length(x)]]
            #counts_per_bam[[group_name]][[sampleID]][["sigi"]]

            #QC filter per BAM
            sigi.filtered1 <- filter_genomicRegions(sigi.ranges, Regions=RegionsToFilter, type="filterOut_overlapping_sites")
            snp.ranges <- sigi.filtered1[[length(sigi.filtered1)]]
            sigi.filtered2 <- filter_genomicRegions(snp.ranges, Regions=RegionsToKeep, type="keep_overlapping_sites")

        #merge results in one list of snp.range that pass filters in each step
        if (is.null(RegionsToFilter))  {res <- counts_per_bam[[group_name]][[sampleID]]}
        if (!is.null(RegionsToFilter) & is.null(RegionsToKeep)) {  res <- c(counts_per_bam[[group_name]][[sampleID]], sigi.filtered1) }
        if (is.null(RegionsToFilter)  & !is.null(RegionsToKeep)) { res <- c(counts_per_bam[[group_name]][[sampleID]], sigi.filtered2) }
        if (!is.null(RegionsToFilter) & !is.null(RegionsToKeep)) { res <- c(counts_per_bam[[group_name]][[sampleID]], sigi.filtered1, sigi.filtered2) }
        res_per_bam[[group_name]][[sampleID]] <- res

        #set progress bar
        if (verbose) {setTxtProgressBar(pb, i)}
        i=i+1
    }
  }
  if (verbose) {close(pb)}
  return(res_per_bam)
}
