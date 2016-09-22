#BaalChIP: filtering functions for QC
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


getdframe <- function(granges) {
    #I had to do this function because as.data.frame is not working for some obscure reason and I am too tired to figure it out
    dframe <- data.frame("seqnames"=as.character(seqnames(granges)),
                      "start"=start(granges),
                      "end"=end(granges),
                      "ID"=values(granges)[["ID"]],
                      "REF"=values(granges)[["REF"]],
                      "ALT"=values(granges)[["ALT"]],
                      "REF.counts"=values(granges)[["REF.counts"]],
                      "ALT.counts"=values(granges)[["ALT.counts"]],
                      "Total.counts"=values(granges)[["Total.counts"]],
                      "Foreign.counts"=values(granges)[["Foreign.counts"]],
                      "AR.counts"=values(granges)[["AR"]],stringsAsFactors = FALSE)
    rownames(dframe) <- names(granges)
    return(dframe)
}


get_mergedcounts <- function(celldata, metadata, includeForeign=FALSE){

    #dataframe
    celldata <- lapply(celldata, getdframe)

    #get the data that will be grouped
    metadata <- split(metadata, metadata$target)

    groupCounts <- lapply(metadata, function (targetGroup) {
        #order by replicate: 1,2,3,...
        targetGroup <- targetGroup[order(targetGroup$replicate_number),,drop=FALSE]

        #Get the counts tables to merge
        data2group <- celldata[targetGroup$sampleID]
        if (!includeForeign) {data2group <- lapply(data2group, subset, select = c("ID","REF.counts","ALT.counts"), drop=FALSE)}
        if (includeForeign) {data2group <- lapply(data2group, subset, select = c("ID","REF.counts","ALT.counts","Foreign.counts"), drop=FALSE)}

        #group tables
        group <- merge(data2group[[1]],data2group[[2]], by=c("ID"))
        if (length(data2group) >= 3) {
                for (d in 3:length(data2group)) {
                    group <- merge(group,data2group[[d]], by=c("ID"))
                }
        }

        #change colnames
        nrrep <- nrow(targetGroup)
        if (!includeForeign) {colname <- unlist(lapply(1:nrrep, function (x) paste0(c("REF.","ALT."),x)))}
        if (includeForeign)  {colname <- unlist(lapply(1:nrrep, function (x) paste0(c("REF.","ALT.","FOREIGN."),x)))}
        colnames(group) <- c("ID", colname)

        #return
        return(group)
    })

    #get unique ids
    ids <- unique(unlist(lapply(groupCounts, "[[", "ID")))

    #add complete data table to "res"
    res <- data.frame("ID"=ids,stringsAsFactors=FALSE)

    for (target_name in names(groupCounts)){
        group <- groupCounts[[target_name]]
        if (nrow(group) > 0) {
            group <- data.frame("score"=1, group, stringsAsFactors=FALSE)
            a <- group[match(ids, group$ID),,drop=FALSE]
            a$ID <- ids
            rownames(a) <- NULL
            a$score[is.na(a$score)] <- 0
            a <- a[,- which(colnames(a)=="ID"),drop=FALSE]
            colnames(a) <- paste0(target_name,".",colnames(a))
            res <- cbind(res, a)
        }
    }

    #debug
    #Are there any columns which are always zero?
    nrsamples <- sum(grepl("score", colnames(res)))
    scores <- res[,grepl("score", colnames(res)),drop=FALSE]
    if (nrsamples == 1) {scores = data.frame(scores, stringsAsFactors=FALSE)}
    if (any(rowSums(scores) == 0)) {stop("Something went wrong! There are rows in the merged table for which no TF overlaps")}

    #return
    return(res)
}


applyMergeResults <- function(samples, res_per_bam, includeForeign=FALSE) {
    cells <- names(res_per_bam)
    res_merged <- lapply(cells, function(x) {list()})
    names(res_merged) <- cells

    for (cellname in cells) {
        lastset <- lapply(res_per_bam[[cellname]], function(x) {return(x[[length(x)]])})
        m1 <- get_mergedcounts(celldata=lastset, metadata=samples[samples$group_name==cellname,,drop=FALSE], includeForeign=includeForeign)
        if (nrow(m1) == 0) {m1 <- data.frame()}
        res_merged[[cellname]][["replicates_merged"]] <- m1
        #messsage(paste0("data frame contains ", nrow(m1), " variants\n"))
    }
    return(res_merged)
}

filter_1allele <- function(mergedcounts) {
    REFsums <- rowSums(mergedcounts[,grepl("REF", colnames(mergedcounts)),drop=FALSE], na.rm=TRUE)
    ALTsums <- rowSums(mergedcounts[,grepl("ALT", colnames(mergedcounts)),drop=FALSE], na.rm=TRUE)
    mergedcounts[(REFsums > 0 & ALTsums > 0),,drop=FALSE]
}

applyFilter1allele <- function(res_merged) {

  for (cellname in names(res_merged)) {
        lastset <- res_merged[[cellname]][[length(res_merged[[cellname]])]]
        m2 <- filter_1allele(lastset)
        res_merged[[cellname]][["Only1Allele"]] <- m2
        #message(paste0(cellname, ": data frame contains ", nrow(m2), " variants (pass filter)\n"))
  }
  return(res_merged)
}







