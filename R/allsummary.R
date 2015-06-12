#BaalChIP: functions for retrieving summary information from BaalChIP class
#Ines de Santiago and Wei Liu (2015)


get_filterstats <- function (sampledata, mergeddata=NULL) {
    
    filtnames1 <- names(sampledata[[1]])
    
    #total number of snp ids per filter (combine across all samples)
    allcounts <- sapply(filtnames1, function(filtname) {
                                                a <- lapply(sampledata, "[[",filtname)
                                                length(unique(unlist(sapply(a,names), use.names=F)))
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
	suppressPackageStartupMessages(require(reshape2))
    suppressPackageStartupMessages(require(doBy))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total))
    filtering_stats$cellname <- rownames(filtering_stats)
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_stats <- melt(filtering_stats, id="cellname")
    data2plot_perc <-  melt(filtering_perc, id="cellname")
    means     <- data.frame(summaryBy(value ~ variable , data2plot_stats, FUN=c(mean)))
    meansPERC <- data.frame(summaryBy(value ~ variable , data2plot_perc, FUN=c(mean)))
    return(cbind(means, "perc"=meansPERC[,2]))
}    
    
applyFilteringStats <- function(res_per_bam, res_merged) {
	if (is.null(res_per_bam)) {return(list())}
    filtering_summary_perSample <- t(do.call("cbind",lapply(res_per_bam, sapply, sapply, length)))
    filtering_stats <- lapply(names(res_per_bam), function(cellname) {
    		sampledata=res_per_bam[[cellname]]
    		if (!is.null(res_merged)) {mergeddata=res_merged[[cellname]]}else{mergeddata=NULL}
    		get_filterstats(sampledata, mergeddata)})
    if (all(sapply(filtering_stats,length) == 0)) {return(list())}
    names(filtering_stats) <- names(res_per_bam)
    filtering_stats <- data.frame(do.call("rbind",filtering_stats))
    average_stats <- get_average_stats(filtering_stats)
    #return(list("filtering_summary_perSample"=filtering_summary_perSample, "filtering_stats"=filtering_stats, "average_stats"=average_stats))
	return(list("filtering_stats"=filtering_stats, "average_stats"=average_stats))
} 

summary_QC <- function(object) {
	if (length(object@mergedCounts) == 0) {res_merged=NULL}else{res_merged =object@mergedCounts} 
	if (length(object@alleleCounts) == 0) {res_per_bam=NULL}else{res_per_bam =object@alleleCounts} 
	
	filtering_stats <- applyFilteringStats(res_per_bam, res_merged)
	
	if(length(filtering_stats)==0) {return(NULL)}
	return(filtering_stats)
}

summary_ASB <- function(object) {
	
	if (length(object@ASB) == 0) {return(NULL)}else{asb <- object@ASB}
	
	asb_stats <- do.call("rbind", lapply(asb, function(x) {c("Ref"=sum(x$Bayes_sig_A),"Alt"=sum(x$Bayes_sig_B))}))
    asb_stats <- cbind(asb_stats, "Total"=rowSums(asb_stats))
    return(asb_stats)
	
}


 
#setGeneric(name="simulationStats",
#                       def=function(.Object,...)
#                       {
#                               standardGeneric("simulationStats")
#                       }
#                       )
#
#setMethod(
#  "simulationStats",
#  "BaalChIP",
#  function(.Object) {
#    FiltIntBias = .Object@param[["QCparam"]][["FiltIntBias"]]
#    if (FiltIntBias) {
#    		query <- Object@QCstats[["simulation_stats"]]
#    		return(query)
#    }else{print("FiltIntBias=FALSE"); return(NULL)}
#  }
#)
