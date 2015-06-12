#BaalChIP: estimation of reference bias funtions
#Ines de Santiago and Wei Liu (2015)

#Refbias functions
getstats <- function(x){
  A <- sum(x$REF.counts, na.rm=T)
  B <- sum(x$ALT.counts, na.rm=T)
  N <- nrow(x)
  return(c("A"=A,"B"=B,"N"=N))
}

#Dont filter if coverage is low due to duplicate removal
filtertable <- function(counts, remove_percentile=0.25) {
    if (remove_percentile>0) {
      #down-sampling reads of sites in the coverage percentile 
      top <- quantile(c(counts$Total.counts), (1-remove_percentile))
      top <- round(top)
      keep <- (counts$Total.counts <= top)
      ignored <- nrow(counts) - sum(keep)
      counts <- counts[ keep,]
      print(paste("ignoring sites with total read number > ",top,": --> ",ignored, "sites ignored during RMbias estimation"))
      return(counts)
    }
    if (remove_percentile==0) {
      print(paste0("keeping all sites"))
      return(counts)
    }
}

pooldata <- function(merged)
    {
    REF.counts <- rowSums(merged[,which(grepl("REF", colnames(merged)))],na.rm=T)
    ALT.counts <- rowSums(merged[,which(grepl("ALT", colnames(merged)))],na.rm=T)
    Total.counts <- REF.counts + ALT.counts
    AR <- REF.counts / Total.counts
    res <- data.frame("ID"=merged$ID, REF.counts, ALT.counts, Total.counts, AR, stringsAsFactors=F)
    return(res)
    }

get_ARestimate <- function(counts, min_n=200) {

  #down-sampling reads of sites in the coverage percentile
  #remove_percentile=0.25
  #counts <- filtertable(counts, remove_percentile=remove_percentile)

  #get global estimation 
  globalEst <- getstats(counts)
  globalEst <- globalEst["A"] / (globalEst["A"] + globalEst["B"])

  #Get estimation per SNP combination
  datasplit <- split(counts, counts$GT)
  datasplit <- lapply(datasplit, getstats)

  #Get AR estimate
  datasplit <- do.call("rbind", datasplit)
  N <- datasplit[,"N"]
  RES <- datasplit[,"A"] / (datasplit[,"A"] + datasplit[,"B"])

  #Order results
  gtord <- c('AG','GA','TC','CT','AC','CA','TG','GT','AT','TA','CG','GC')
  RES <- RES[gtord]
  N <- N[gtord]
  N[is.na(N)] <- 0
  names(RES) <- names(N) <- gtord

  #Replace by the global estimate if N<min_n
  toreplace <- names(RES[N < min_n])
  RES[N < min_n] <- globalEst

  return(RES)
}

estimateRefBias <- function(assayed, GTtable, min_n=200) {

  #Pool data
  counts <- pooldata(assayed)

  #merge tables
  counts <- merge(GTtable, counts, by="ID")
  if(nrow(counts) != nrow(assayed)) {stop("an error occured")}

  #Get genotype for snps in pooled data
  GT <- paste0(counts$REF, counts$ALT)
  counts$GT <- GT

  #Get only SNPs with RAF_tr
  #if(!is.null(RAF_tr)) {
  #  keep <- counts$RAF < RAF_tr[2] & counts$RAF > RAF_tr[1]
  #  ignored <- nrow(counts) - sum(keep)
  #  print(paste("keeping only sites with RAF within [",RAF_tr[1],",",RAF_tr[2],"] --> ",ignored, "sites ignored during RMbias estimation"))
  #  counts <- counts[keep,]
  #}

  #estimate ref bias
  #ARestimate <- get_ARestimate(counts, remove_percentile=remove_percentile, min_n=min_n)
  ARestimate <- get_ARestimate(counts, min_n=min_n)
  print("ARestimate:")
  print(ARestimate)
  ARestimate
}

getbiasTable <- function(assayed, GTtable, ARestimate){

  if (!("RAF" %in% colnames(GTtable))) {GTtable$RAF <- 0.5}
  biastable <- GTtable[GTtable$ID %in% assayed$ID,c("ID","RAF")]
  
  if (!is.null(ARestimate)) {
    GT <- paste(GTtable$REF, GTtable$ALT, sep="")
    biastable$RMbias <- ARestimate[match(GT, names(ARestimate))]
  }else{
    biastable$RMbias <- 0.5
  }
    
  biastable <- biastable[,c("ID","RMbias", "RAF")]
  
  #check
  if (nrow(assayed) != nrow(biastable)) {stop('some error occured')}
  
  #order
  m <- merge(assayed, biastable, by="ID")
  assayed <- m[,1:ncol(assayed)]
  biastable <- m[,c(1, (ncol(assayed)+1):ncol(m))]
  if (any(assayed$ID != biastable$ID)) {stop('some error occured')}

  return(list(assayed, biastable))
}


getbiasparam <- function(biastable){
    biasparam = c("RAF"=TRUE,"RMbias"=TRUE)
    if (all(biastable$RMbias == 0.5)) {biasparam["RMbias"] <- FALSE}
    if (all(biastable$RAF == 0.5)) {biasparam["RAF"] <- FALSE}
    return (biasparam)
}

#get_assayedVar 
get_Vartable <- function(assayedVar, hets) {
    Vartable <- list()
    for (ID in names(hets)){
        assayed <- assayedVar[[ID]]
        snps <- read.delim(hets[[ID]], stringsAsFactors=F, head=T) 
        snps <- snps[snps$ID %in% assayed$ID,]
        Vartable[[ID]] <- snps
    }
    return(Vartable)
}
