#BaalChIP: filtering functions for QC
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


filter_regions <- function(snp.ranges, region.ranges, type=c("filterOut_overlapping_sites","keep_overlapping_sites")) {
    suppressPackageStartupMessages(require(GenomicRanges))
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


applyFiltersPerBam <- function(counts_per_bam, RegionsToFilter, RegionsToKeep) {

    res_per_bam <- list()
    
    N <- sum(sapply(counts_per_bam, length))
    
    cat("-applying filters per BAM\n")
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    i=1
    
    for (cellname in names(counts_per_bam)) {
    	for (sampleID in names(counts_per_bam[[cellname]])) {

			#GRanges of variants		
    		sigi.ranges <- counts_per_bam[[cellname]][[sampleID]][["sigi"]]
    		
    		#QC filter per BAM
    		sigi.filtered1 <- filter_genomicRegions(sigi.ranges, Regions=RegionsToFilter, type="filterOut_overlapping_sites")
    		snp.ranges <- sigi.filtered1[[length(sigi.filtered1)]]
    		sigi.filtered2 <- filter_genomicRegions(snp.ranges, Regions=RegionsToKeep, type="keep_overlapping_sites")
    
            #merge results in one list of snp.range that pass filters in each step
            if (is.null(RegionsToFilter))  {res <- list("sigi"=sigi.ranges)}
            if (!is.null(RegionsToFilter) & is.null(RegionsToKeep)) {  res <- c(list("sigi"=sigi.ranges), sigi.filtered1) }
            if (is.null(RegionsToFilter)  & !is.null(RegionsToKeep)) { res <- c(list("sigi"=sigi.ranges), sigi.filtered2) }
            if (!is.null(RegionsToFilter) & !is.null(RegionsToKeep)) { res <- c(list("sigi"=sigi.ranges), sigi.filtered1, sigi.filtered2) }
            res_per_bam[[cellname]][[sampleID]] <- res    
            
            #set progress bar
        	setTxtProgressBar(pb, i)
        	i=i+1
        }
    }
    close(pb)
    return(res_per_bam)
}