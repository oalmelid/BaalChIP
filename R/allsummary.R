#BaalChIP: functions for retrieving summary information from BaalChIP class
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


get_filterstats <- function (sampledata, mergeddata=NULL) {

    filtnames1 <- names(sampledata[[1]])

    #total number of snp ids per filter (combine across all samples)
    allcounts <- sapply(filtnames1, function(filtname) {
                                                a <- lapply(sampledata, "[[",filtname)
                                                length(unique(unlist(lapply(a,names), use.names=FALSE)))
                                    })
    names(allcounts) <- filtnames1

    #total number of snp ids per filter (for merged counts)
    if (!is.null(mergeddata)) {allcounts <- c(allcounts, sapply(mergeddata, nrow))}

    #number of filtered
    n <- length(allcounts)
    if(n==1) {return(NULL)} #nothing was filtered
    nrfilt <- allcounts[1: (n-1) ] - allcounts[2: (n) ]
    names(nrfilt) <- names(allcounts)[2:(n)]
    nrfilt <- c(nrfilt, "pass"=as.numeric(allcounts[n]))

    #check
    total = allcounts[1]
    if (sum(nrfilt) != total) { stop("error counting filtered SNPs", call.=FALSE) }


    return(nrfilt)
}


get_average_stats <- function(filtering_stats) {
    #suppressPackageStartupMessages(require(reshape2))
    #suppressPackageStartupMessages(require(doBy))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total), stringsAsFactors=FALSE)
    filtering_stats$cellname <- rownames(filtering_stats)
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_stats <- melt(filtering_stats, id="cellname")
    data2plot_perc <-  melt(filtering_perc, id="cellname")
    means <- data.frame(summaryBy(value ~ variable , data2plot_stats, FUN=c(mean)), stringsAsFactors=FALSE)
    meansPERC <- data.frame(summaryBy(value ~ variable , data2plot_perc, FUN=c(mean)), stringsAsFactors=FALSE)
    return(cbind(means, "perc"=meansPERC[,2]))
}

applyFilteringStats <- function(res_per_bam, res_merged) {
    if (is.null(res_per_bam)) {return(list())}
    filtering_summary_perSample <- t(do.call("cbind",lapply(res_per_bam, sapply, sapply, length)))
    filtering_stats <- lapply(names(res_per_bam), function(cellname) {
            sampledata=res_per_bam[[cellname]]
            if (!is.null(res_merged)) {mergeddata=res_merged[[cellname]]}else{mergeddata=NULL}
            if (all(sapply(mergeddata,length) == 0)) {mergeddata <- NULL}
            get_filterstats(sampledata, mergeddata)
            })

    if (all(sapply(filtering_stats,length) == 0)) {return(list())}
    names(filtering_stats) <- names(res_per_bam)
    filtering_stats <- data.frame(do.call("rbind",filtering_stats), stringsAsFactors=FALSE)
    average_stats <- get_average_stats(filtering_stats)
    #return(list("filtering_summary_perSample"=filtering_summary_perSample, "filtering_stats"=filtering_stats, "average_stats"=average_stats))
    return(list("filtering_stats"=filtering_stats, "average_stats"=average_stats))
}

summary_QC <- function(object) {
    mergedCounts <- getBaalSlot(object, "mergedCounts")
    alleleCounts <- getBaalSlot(object, "alleleCounts")
    if (length(mergedCounts) == 0) {res_merged=NULL}else{res_merged = mergedCounts}
    if (length(alleleCounts) == 0) {res_per_bam=NULL}else{res_per_bam = alleleCounts}

    filtering_stats <- applyFilteringStats(res_per_bam, res_merged)

    if(length(filtering_stats)==0) {
        message("-no filters applied yet")
        return(NULL)
    }
    return(filtering_stats)
}

summary_ASB <- function(object) {

    asb <- getBaalSlot(object, "ASB")
    if (length(asb) == 0) {return(NULL)}
    asb_stats <- do.call("rbind", lapply(asb, function(x) {c("Ref"=sum(x$Bayes_sig_A),"Alt"=sum(x$Bayes_sig_B))}))
    asb_stats <- cbind(asb_stats, "Total"=rowSums(asb_stats))
    return(asb_stats)

}


getReport <- function(object, group_name) {

    asb <- getBaalSlot(object, "ASB")
    if (length(asb) == 0) {return(NULL)}

    asb <- asb[[group_name]]
    assayed <- getBaalSlot(object, "assayedVar")[[group_name]]
    pooled <- pooldata(assayed)
    VarTable <- getBaalSlot(object, "VarTable")[[group_name]][,c("ID","CHROM","POS","REF","ALT")]
    RAFtable <- getBaalSlot(object, "biasTable")[[group_name]]
    correctedAR <- data.frame("ID"=as.character(asb$ID),
                    "Bayes_lower"=asb$Bayes_lower,
                    "Bayes_upper"=asb$Bayes_upper,
                    "Corrected.AR"=rowMeans(asb[,c("Bayes_lower","Bayes_upper")]),
                     stringsAsFactors=FALSE)
    baalSig <- as.character(asb$ID[asb[,"Bayes_sig_A"]==1 | asb[,"Bayes_sig_B"] == 1])
    snps <- merge(VarTable, pooled, by="ID")
    snps <- merge(snps, RAFtable, by="ID")
    snps <- merge(snps, correctedAR, by="ID")
    snps$isASB <- snps$ID %in% baalSig
    return(snps)
}


