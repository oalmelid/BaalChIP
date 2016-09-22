#BaalChIP: estimation of reference bias ans ASB funtions
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

pooldata <- function(merged) {
    REF.counts <- rowSums(merged[,which(grepl("REF", colnames(merged)))],na.rm=TRUE)
    ALT.counts <- rowSums(merged[,which(grepl("ALT", colnames(merged)))],na.rm=TRUE)
    Total.counts <- REF.counts + ALT.counts
    AR <- REF.counts / Total.counts
    res <- data.frame("ID"=merged$ID, REF.counts, ALT.counts, Total.counts, AR, stringsAsFactors=FALSE)
    return(res)
    }

estimateRefBias <- function(assayed, GTtable, min_n=200) {

    if (nrow(assayed) == 0) {return(data.frame())} #in case there were no variants left

    #Pool data
    counts <- pooldata(assayed)
    globalEst <- mean(counts$AR)

    #merge tables
    counts <- merge(GTtable, counts, by="ID")
    if(nrow(counts) != nrow(assayed)) {stop("an error occured")}

    #Get genotype for snps in pooled data
    GT <- paste0(counts$REF, counts$ALT)
    counts$GT <- GT

    ##Get AR estimate
    datasplit <- split(counts, counts$GT)
    RES <- lapply(datasplit, function(x) {mean(x$AR)})
    N <- unlist(lapply(datasplit, nrow))

    #Order results
    gtord <- c('AG','GA','TC','CT','AC','CA','TG','GT','AT','TA','CG','GC')
    RES <- RES[gtord]
    N <- N[gtord]
    RES[is.na(names(RES))] <- NA
    N[is.na(N)] <- 0
    names(RES) <- names(N) <- gtord

    #Replace by the global estimate if N<min_n
    toreplace <- names(RES[N < min_n])
    RES[N < min_n] <- globalEst

    unlist(RES)
}

#' @importFrom stats na.omit
getbiasTable <- function(assayed, GTtable, ARestimate){

    if(nrow(assayed)==0) {return(list(data.frame(), data.frame()))}

    biastable <- GTtable[GTtable$ID %in% assayed$ID,c("ID","RAF"), drop=FALSE]

    if (!is.null(ARestimate)) {
    GT <- paste(GTtable$REF, GTtable$ALT, sep="")
    biastable$RMbias <- ARestimate[match(GT, names(ARestimate))]
    }else{
    biastable$RMbias <- 0.5
    }

    biastable <- biastable[,c("ID","RMbias", "RAF"), drop=FALSE]

    #check
    if (nrow(assayed) != nrow(biastable)) {stop('some error occured')}
    if (!(is.numeric(biastable[,"RMbias"]))) {stop('RM is not numeric')}
    if (!(is.numeric(biastable[,"RAF"]))) {stop('RAF is not numeric')}

    #order
    biastable <- data.frame(na.omit(biastable))
    m <- merge(assayed, biastable, by="ID")
    if (nrow(m) == 0) {stop('No variants left after merging counts with biastable')}
    assayed <- m[,1:ncol(assayed), drop=FALSE]
    biastable <- m[,c(1, (ncol(assayed)+1):ncol(m)), drop=FALSE]
    if (any(assayed$ID != biastable$ID)) {stop('some error occured')}

    return(list(assayed, biastable))
}


getbiasparam <- function(biastable){
    if (nrow(biastable) == 0) {
        return(c("RAF"=NA, "RMbias"=NA))
    }

    biasparam = c("RAF"=TRUE,"RMbias"=TRUE)
    if (all(biastable$RMbias == 0.5)) {biasparam["RMbias"] <- FALSE}
    if (all(biastable$RAF == 0.5)) {biasparam["RAF"] <- FALSE}
    return (biasparam)
}


getRAFfromgDNA <- function (bamFiles, snp.ranges, min_base_quality=10, min_mapq=15, verbose=TRUE) {
    if (verbose) {
        message("-computing allele counts per gDNA BAM")
        pb <- txtProgressBar(min = 0, max = length(bamFiles), style = 3)
    }

    AllCounts <- lapply(seq_len(length(bamFiles)), function(i) {
        bamfile <- bamFiles[i]
        acounts <- get_allele_counts(bamfile, snp.ranges, returnRanges=FALSE, min_base_quality=min_base_quality,min_mapq=min_mapq)
        acounts <- acounts[,c("ID","CHROM","POS","REF.counts","ALT.counts"), drop=FALSE]
        if (verbose) {setTxtProgressBar(pb, i)}
        return(acounts)
    })

    names(AllCounts) <- bamFiles
    if (verbose) {close(pb)}

    #delete all entries in AllCounts that are NULL
    if (any(sapply(AllCounts, is.null))) {
        idx <- which(!(sapply(AllCounts, is.null)))
        AllCounts <- AllCounts[idx]
    }

    snpIDs <- names(snp.ranges)
    if (length(AllCounts) > 1) {
        refcounts <- rowSums(sapply(AllCounts, function(x) {x$REF.counts[match(snpIDs, x$ID)]}), na.rm=TRUE)
        altcounts <- rowSums(sapply(AllCounts, function(x) {x$ALT.counts[match(snpIDs, x$ID)]}), na.rm=TRUE)
    }

    if (length(AllCounts) == 1)   {
        x <- AllCounts[[1]]
        refcounts <- x$REF.counts[match(snpIDs, x$ID)]
        altcounts <- x$ALT.counts[match(snpIDs, x$ID)]
    }

    if (length(AllCounts) == 0) {
        return(NULL)
    }

    refcounts [ is.na(refcounts) ] <- 0
    altcounts [ is.na(altcounts) ] <- 0
    totalcounts <- (refcounts + altcounts)
    RAF <- refcounts / totalcounts
    RAF <- data.frame("ID"=snpIDs, "RAF"=RAF, stringsAsFactors=FALSE)
    #RAF[totalcounts <= 5] <- NA
    return(RAF)

}

useRAFfromhets <- function(snps, ID, verbose=TRUE) {
    if (nrow(snps)==0) {return(data.frame())} #in case there were no snps left after filtering

    if (!("RAF" %in% colnames(snps))) {
        #oops... there are no RAF value in hets table...
        warning("no RAF values found for ",ID," will use 0.5 for all variants")
        snps$RAF <- 0.5
    }else{
        #will use the RAF values in this case:
        if (verbose) {message("will use RAF values for ", ID," group from hets file")}
    }
        return(snps)
}

useRAFfromgDNA <- function(gDNAbams, snps, ID, min_base_quality=10, min_mapq=15, verbose=TRUE) {

    if (nrow(snps)==0) {return(data.frame())} #in case there were no snps left after filtering

    if (verbose) {message("-calculating RAF from gDNA for group ",ID)}
    bamFiles <- gDNAbams
    snp.ranges <- get_snp_ranges(snps)
    RAF <- getRAFfromgDNA(bamFiles, snp.ranges, min_base_quality=min_base_quality, min_mapq=min_mapq, verbose=verbose)
    if (!is.null(RAF)) {snps <- merge(snps, RAF, by="ID")}
    if (is.null(RAF)) {stop("Did not find any reads from gDNA BAMs overlapping SNPs for group",ID," Cannot proceed\n")}
    #message("will use RAF values estimated from all gDNA bam files for ",ID, " group:", gDNA[[ID]])
    return(snps)
}

get_Vartable <- function(assayedVar, hets, gDNA=list(), min_base_quality=10, min_mapq=15, RAFcorrection=TRUE, verbose=TRUE) {

    if (length(gDNA)==0) {gDNA <- NULL}

    Vartable <- lapply(names(hets), function (ID) {
        #for (ID in names(hets)){
        snps <- read.delim(hets[[ID]], stringsAsFactors=FALSE, header=TRUE)
        assayed <- assayedVar[[ID]]
        snps <- snps[snps$ID %in% assayed$ID,,drop=FALSE]
        gDNAbams <- gDNA[[ID]]


        #RAF correction is TRUE, and gDNA is null --> go directly to RAF
        if (is.null(gDNAbams) & RAFcorrection) { snps <- useRAFfromhets(snps, ID, verbose=verbose) }

        #RAF correction is TRUE and gDNA is not null
        if (!is.null(gDNAbams) & RAFcorrection) {

            if ("RAF" %in% colnames(snps)) {
                #There are both gDNA and RAF in hets tables.. will use the RAF instead!
                warning("both gDNA and hets file found for group ", ID, ". Will use RAF from hets!")
                snps <- useRAFfromhets(snps, ID, verbose=verbose)
            }else{
                snps <- useRAFfromgDNA(gDNAbams, snps, ID, min_base_quality=min_base_quality, min_mapq=min_mapq, verbose=verbose)

            }


        }

    #Vartable[[ID]] <- snps
    return(snps)
    })
    names(Vartable) <- names(hets)

    return(Vartable)
}


addVarTable <- function(object, RAFcorrection=TRUE, verbose=TRUE) {
    assayedVar <- getBaalSlot(object, "assayedVar")
    hets <- getBaalSlot(object, "hets")

    ##-----Read hets table and update "RAF" values accordingly
    QCparam <- getBaalSlot(object, "param")$QCparam
    gDNA=getBaalSlot(object, "gDNA")
    min_base_quality <- QCparam[["min_base_quality"]] #will not be used if gDNA is null or RAFcorrection==FALSE
    min_mapq <- QCparam[["min_mapq"]] #will not be used if gDNA is null or RAFcorrection==FALSE

    VarTable <- get_Vartable(assayedVar, hets, gDNA, min_base_quality, min_mapq, RAFcorrection,verbose)
    return(VarTable)
}
