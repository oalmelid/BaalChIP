#BaalChIP: estimation of reference bias ans ASB funtions
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

pooldata <- function(merged)
    {
    REF.counts <- rowSums(merged[,which(grepl("REF", colnames(merged)))],na.rm=T)
    ALT.counts <- rowSums(merged[,which(grepl("ALT", colnames(merged)))],na.rm=T)
    Total.counts <- REF.counts + ALT.counts
    AR <- REF.counts / Total.counts
    res <- data.frame("ID"=merged$ID, REF.counts, ALT.counts, Total.counts, AR, stringsAsFactors=F)
    return(res)
    }

estimateRefBias <- function(assayed, GTtable, min_n=200) {
  
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
