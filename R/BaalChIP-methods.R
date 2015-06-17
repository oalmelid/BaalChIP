#BaalChIP: all methods
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

##initialization method
setMethod("initialize",
          "BaalChIP",
          function(.Object, samplesheet = NULL, hets=NULL) {
					
            ##-----check arguments
            if(missing(samplesheet))stop("NOTE: 'samplesheet' is missing!")    
            samples <- BaalChIP.checks(name="samplesheet",samplesheet)
            BaalChIP.checks(name="hets",hets)
            
            ##-----check matching cellnames
            cells1 <- unique(.Object@samples$cell_name)
            checkmatchingnames(names(.Object@hets), cells1)
    
            ##-----initialization
            .Object@samples <- samples
            .Object@hets <- hets
            .Object@param <- list()
            			
            .Object
          }
)

setMethod("show", "BaalChIP",
    function(object){
    cat(" Type :", class(object), "\n")
    QCstats <- summaryQC(object)[["filtering_stats"]]
    asb <- object@ASB
    samples <- object@samples
    cat(" Samples                 :  ", nrow(samples), "\n")
    cat(" Experiments             :  ", unique(samples[,"cell_name"]), "\n")
    cat(" Filtering and QC        :  ", ifelse(is.null(QCstats), "None", paste(ncol(QCstats), "filters applied")), "\n")
    cat(" Run allele-specific     :  ", ifelse(length(asb)==0, "None", "Yes: run BaalChIP.get(object, 'ASB')"), "\n")
    cat("\n")
})

#' Generates allele-specific read count data
#' @name alleleCounts
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @description Generates allele-specific read count data from each BAM ChIP-seq dataset for each variant.
#' @param .Object An object of the \code{\link{BaalChIP}} class 
#' @param min_base_quality A numeric value indicating the minumum read base quality below which the base is ignored when summarizing pileup information (default 10)
#' @param min_mapq A numeric value indicating the minumum mapping quality (MAPQ) below which the entire read is ignored (default 15)
#' @note BaalChIP computes allelic counts at each variant position with Rsamtools pileup function. The algorithm follows pileup::Rsamtools by automatically excluding reads flagged as unmapped, secondary, duplicate, or not passing quality controls.
#' @details Utilizes the information whithin the \code{samples} slot of a BaalChIP object. Will primarly find all variants overlapping peaks. Then, for each variant, computes the number of reads carrying the reference (REF) and alternative (ALT) alleles.
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{alleleCounts} containing a list of GRanges objects.
#' @seealso \code{\link{BaalChIP.get}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'
#'#retrieve alleleCounts:
#'counts <- BaalChIP.get(res, "alleleCountsPerBam")
#'
#'#alleleCounts are grouped by bam_name and group_name:
#'names(counts)
#'names(counts[["MCF7"]])
#'
#'#check out the result for one of the bam files:
#'counts[["MCF7"]][[1]]
#' @export 
setMethod(
  f="alleleCounts",
  signature="BaalChIP",
  function(.Object, 
		   min_base_quality=10, 
		   min_mapq=15){
								
  ##-----check input arguments
  samples <- .Object@samples
  hets <- .Object@hets
  BaalChIP.checks(name="min_base_quality", min_base_quality)
  BaalChIP.checks(name="min_mapq", min_mapq)
   
  ##-----assign parameters
  .Object@param$QCparam <- list(min_base_quality=min_base_quality, min_mapq=min_mapq)
		   
  ##-----compute allele counts
  res_per_bam  <- applyAlleleCountsPerBam(samples, hets, min_base_quality, min_mapq)
  .Object@alleleCounts <- res_per_bam
  
  return(.Object)
  }
)
    
#' Removes variants that may be problematic for identification of allele-specific events
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name QCfilter
#' @description Quality control step for removing variants that may be problematic for identification of allele-specific events.
#' @param .Object An object of the \code{\link{BaalChIP}} class 
#' @param RegionsToFilter a named list of GRanges objects. Variants overlapping these regions will be removed
#' @param RegionsToKeep a named list of GRanges objects. Works in an oposite way to 'RegionstoFilter', variants NOT overlapping these regions will be removed
#' @return An updataded \code{\link{BaalChIP}} object with the slot \code{alleleCounts} containing a list of GRanges objects that pass filters.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data(blacklist_hg19)
#'data(pickrell2011cov1_hg19)
#'data(UniqueMappability50bp_hg19)
#'res <- QCfilter(res, 
#'                RegionsToFilter=list("blacklist"=blacklist_hg19, "highcoverage"=pickrell2011cov1_hg19), 
#'                RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
#'
#'#retrieve alleleCounts:
#'counts <- BaalChIP.get(res, "alleleCountsPerBam")
#'
#'#alleleCounts are grouped by bam_name and group_name:
#'names(counts)
#'names(counts[["MCF7"]])
#'
#'#check out the result for one of the bam files:
#'counts[["MCF7"]][[1]]
#' @export 
setMethod(
  f="QCfilter",
  signature="BaalChIP",
  function(.Object, 
           RegionsToFilter=list("blacklist"=blacklist_hg19, "highcoverage"=pickrell2011cov1_hg19),
           RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19)){
								
  ##-----check input arguments
  RegionsToFilter <- BaalChIP.checks(name="RegionsToFilter", RegionsToFilter)
  RegionsToKeep <- BaalChIP.checks(name="RegionsToKeep", RegionsToKeep)
    
  ##-----assign parameters
  counts_per_bam <- .Object@alleleCounts
  .Object@param$QCparam$RegionsToFilter <- names(RegionsToFilter)
  .Object@param$QCparam$RegionsToKeep <- names(RegionsToKeep)
		   
  ##-----apply filters
  res_per_bam <- applyFiltersPerBam(counts_per_bam, RegionsToFilter, RegionsToKeep)
  .Object@alleleCounts <- res_per_bam
  
  return(.Object)
  }
)

#setMethod(
#  f="filterIntbias",
#  signature="BaalChIP",
#  function(.Object, simul_output){
#								
#  ##-----check input arguments
#  BaalChIP.checks(name="simul_output", simul_output)
#    
#  ##-----assign parameters
#  res_per_bam <- .Object@alleleCounts 
#  samples <- .Object@samples
#  .Object@param$QCparam$FiltIntBias <- TRUE
#  .Object@param$QCparam$simul_output <- simul_output
#	
#  ##-----read length per sample
#  samples$readlen <- applyReadlenPerBam(samples, res_per_bam)
#   	   
#  ##---- run simulations
#  if (testingScript) {
#    	#save so far before running simulations
#    	tmpname <- file.path(dirname(samplesheet),paste0(tmpfile_prefix,"_step1.Rda"))
#    	save(samples, res_per_bam, file = tmpname)
#  }
#    
#  simul_output <- file.path(simul_output,tmpfile_prefix)
#  simures <- applySim(samples, res_per_bam, simul_output=simul_output, simulation_script=simulation_script, testingScript=testingScript)
#  
#  ##----apply filter
#  res_per_bam <- applyIntBiasFilterPerBam(samples, res_per_bam, simcounts=simures[["simcounts"]])
#  
#  # save so far
#  if (testingScript) {
#    	tmpname <- file.path(dirname(samplesheet),paste0(tmpfile_prefix,"_step2.Rda"))
#    	save(samples, res_per_bam, simures, file = tmpname)
#  }
#   
#  
#  .Object@alleleCounts <- res_per_bam
#  return(.Object)
#  }
#)


#' Merges allele-specific read count data per group
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name mergePerGroup
#' @description Merges all ChIP-seq datasets within a group of samples creating a data.frame that contains allele-specific read count data for all variants that need to be analysed. 
#' @param .Object An object of the \code{\link{BaalChIP}} class 
#' @details if QCfilter has been applied, will use the most up-to-date variant set available for each individual BAM file (after QC). Missing values are allowed for heterozygous variants that are not available (e.g. do not pass filter for a particular ChIP-seq dataset).
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{mergedCounts} containing a data.frame of merged samples per group.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data(blacklist_hg19)
#'data(pickrell2011cov1_hg19)
#'data(UniqueMappability50bp_hg19)
#'res <- QCfilter(res, 
#'                RegionsToFilter=list("blacklist"=blacklist_hg19, "highcoverage"=pickrell2011cov1_hg19), 
#'                RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
#'
#'res <- mergePerGroup(res)
#'
#'#retrieve mergedCounts:
#'counts <- BaalChIP.get(res, "mergedCounts")
#'
#'#mergedCounts are grouped by group_name:
#'names(counts)
#'sapply(counts, dim)
#'
#'#check out the result for one of the groups:
#'head(counts[[1]])
#' @export 
setMethod(
  f="mergePerGroup",
  signature="BaalChIP",
  function(.Object){
  							
  ##-----assign parameters
  res_per_bam <- .Object@alleleCounts
  samples <- .Object@samples
  
  ##-----apply filters
  res_merged <- applyMergeResults(samples, res_per_bam)
  .Object@mergedCounts <- res_merged
  
  return(.Object)
  }
)
#' Filters out variants with only 1 observed allele
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name filter1allele
#' @description Filters the data frame available within a \code{\link{BaalChIP}} object (slot \code{mergedCounts}). This filter ignores variants for which only one allele is observed after pooling ChIP-seq reads from all datasets. 
#' @param .Object An object of the \code{\link{BaalChIP}} class 
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{mergedCounts} containing a data.frame of merged samples per group with variants that pass the filter.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data(blacklist_hg19)
#'data(pickrell2011cov1_hg19)
#'data(UniqueMappability50bp_hg19)
#'res <- QCfilter(res, 
#'                RegionsToFilter=list("blacklist"=blacklist_hg19, "highcoverage"=pickrell2011cov1_hg19), 
#'                RegionsToKeep=list("UniqueMappability"=UniqueMappability50bp_hg19))
#'
#'res <- mergePerGroup(res)
#'res <- filter1allele(res)
#'
#'#retrieve mergedCounts:
#'counts <- BaalChIP.get(res, "mergedCounts")
#'
#'#mergedCounts are grouped by group_name:
#'names(counts)
#'sapply(counts, dim)
#'
#'#check out the result for one of the groups:
#'head(counts[[1]])
#' @export 
setMethod(
   f="filter1allele",
   signature="BaalChIP",
   function(.Object){
  							
   ##-----assign parameters
   res_merged <- .Object@mergedCounts
  
   ##-----apply filters
   res_merged <- applyFilter1allele(res_merged)
   .Object@mergedCounts <- res_merged
  
   return(.Object)
   }
)



#' BaalChIP pipeline - allele counts and QC
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name BaalChIP.QC
#' @description BaalChIP.QC is a wrapper convenience function, to compute allele counts and perform quality controls in one step. This function will use the package's defaults.
#' @param .Object An object of the \code{\link{BaalChIP}} class
#' @details This function is a wrapper of the following functions: \code{\link{alleleCounts}}, \code{\link{QCfilter}}, \code{\link{mergePerGroup}}, \code{\link{filter1allele}}
#' @return An object of the \code{\link{BaalChIP}} class.
#' @seealso \code{\link{summaryQC}}, \code{\link{plotQC}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- BaalChIP.QC(res)
#'
#'#summary of the QC step
#'summaryQC(res)
#' @export 

setMethod(
  f="BaalChIP.QC",
  signature="BaalChIP",
  function(.Object){
		#simul_output=NULL)
								
    ##----- compute allele counts
    .Object <- alleleCounts(.Object, min_base_quality=10, min_mapq=15)
	
    ##-----run QC step
    data(blacklist_hg19)
    data(pickrell2011cov1_hg19)
    data(UniqueMappability50bp_hg19)
    .Object <- QCfilter(.Object, RegionsToFilter=c("blacklist_hg19","pickrell2011cov1_hg19"), 
                        RegionsToKeep=c("UniqueMappability50bp_hg19"))
    
    #if (FiltIntBias == FALSE & !is.null(simul_output)) {
	  #		warning (paste("will not use 'simul_output' because FiltIntBias is FALSE"))
	  #
	
    ##-----merge replicates
    cat("-merging replicated samples\n")
    .Object <- mergePerGroup(.Object)
    
    ##-----filter 'Only1Allele'
    cat("-filtering out SNPs with only 1 observed allele\n")
    .Object <- filter1allele(.Object)
    
    cat("-QC complete!\n\n")
    .Object
  }
)

#' Identifies allele-specific binding events
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name getASB
#' @description getASB identifies allele-specific binding events using a bayesian framework.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param Iter Maximun number of iterations (default 5000)
#' @param conf_level Confidence interval in the estimated allelic ratio (default 0.95)
#' @param RMcorrection Logical value indicating if reference mapping (RM) bias should be applied (default TRUE). If FALSE will not correct for reference allele mapping bias. If TRUE will estimate the RM bias from the overall reference allele proportion 
#' @param RAFcorrecion Logical value indicating if relative allele frequency (RAF) bias correction should be applyed (default TRUE). If TRUE will read RAF values for each variant from \code{hets} files (RAF column name). If FALSE will not correct for relative allele frequency bias
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{ASB} containing variants identified as allele-specific.
#' @seealso \code{\link{summaryASB}}, \code{\link{BaalChIP.report}}
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)
#'res <- BaalChIP.QC(res)
#'res <- getASB(res)
#'
#'#summary - number of significant ASB variants
#'summaryASB(res)
#'
#'#report result
#'res <- BaalChIP.report(res)
#' @export 
setMethod(
  "getASB",
  "BaalChIP",
  function(.Object, Iter=5000, conf_level=0.95, RMcorrection = TRUE, RAFcorrection=TRUE){
								
    ##-----check input arguments
    #
    BaalChIP.checks(name="Iter", Iter)
    BaalChIP.checks(name="conf_level", conf_level)
    BaalChIP.checks(name="RMcorrection", RMcorrection)
    BaalChIP.checks(name="RAFcorrection", RAFcorrection)  
    
	##-----updade object with "assayedVar" and "GTtable"
	assayedVar <- BaalChIP.get(.Object, "mergedCounts") #last filtered mergedCounts
	hets=.Object@hets
    .Object@assayedVar <- assayedVar
    VarTable <- get_Vartable(assayedVar, hets)
    .Object@VarTable <- VarTable
    
    ##-----check matching cellnames
    cells1 <- unique(.Object@samples$cell_name)
    checkmatchingnames(names(.Object@assayedVar), cells1)
    checkmatchingnames(names(.Object@VarTable), cells1)
    
    ##-----run
    Expnames <- names(assayedVar)
    results <- list()
    biasTable <- list()
    applyedCorrection <- list()
    
    for (ID in Expnames) {
      print (paste("... running for:", ID))
    	assayed <- assayedVar[[ID]]
    	GTtable <- VarTable[[ID]]
    	
    	#if no RAF replace any RAF values by 0.5
    	if (!RAFcorrection) {GTtable$RAF <- 0.5}
    	
    	#get bias table (variants ordered equally between counts table and biastable)
        if (RMcorrection) {ARestimate <- estimateRefBias(assayed, GTtable, min_n=200)}else{ARestimate=NULL}
        result <- getbiasTable(assayed, GTtable, ARestimate) 
    	counts <- result[[1]]
    	biastable <- result[[2]]
    	biasparam <- getbiasparam(biastable)
    	
    	#run bayes
    	Bayes_report <- runBayes(counts=counts, bias=biastable, Iter=Iter, conf_level=conf_level) 
    	
    	#append results
    	results[[ID]] <- Bayes_report
    	biasTable[[ID]] <- biastable	
    	applyedCorrection[[ID]] <- biasparam
    }
    

    ##-----assign parameters
    applyedCorrection <- t(do.call("rbind",applyedCorrection))
	.Object@param$ASBparam <- list(Iter=Iter, conf_level=conf_level, applyedCorrection=applyedCorrection)
	
    ##-----updade status and return
    .Object@biasTable <- biasTable
    .Object@ASB <- results
    
    cat("-ASB identification complete!\n\n")
    .Object
  }
)			

#' Get slots from a BaalChIP object
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name BaalChIP.get
#' @description Get information from individual slots in a BaalChIP object.
#' @param .Object An object of the \code{\link{BaalChIP}} class
#' @param what a single character value specifying which information should be retrieved. Options: 'samples', 'param', 'alleleCountsPerBam', 'mergedCounts'
#' @return The slot content from an object of the \code{\link{BaalChIP}} class.
#' @examples
#'data(BaalObject)
#'
#'#samples data spreadsheet and hets:
#'BaalChIP.get(BaalObject,"samples")
#'
#'#parameters used within run:
#'BaalChIP.get(BaalObject,"param")
#'
#'#retrieve a GRanges list with allele-specific read counts per BAM file:
#'counts <- BaalChIP.get(BaalObject,"alleleCountsPerBam")
#'counts[["MCF7"]][[1]]
#'
#'#retrieve a data.frame with allele-specific read counts per group:
#'counts <- BaalChIP.get(BaalObject,"mergedCounts")
#'head(counts[[1]])
#' @export 
setMethod(
  "BaalChIP.get",
  "BaalChIP",
  function(.Object, what=c("samples","param","alleleCountsPerBam","mergedCounts")) {
    ##-----check input arguments
    BaalChIP.checks(name="get.what",param=what)
    ##-----get query
    query<-NULL
    if(what=="samples"){
      query<- list("samples"=.Object@samples,"hets"=.Object@hets)
    } else if(what=="param"){
      query<-.Object@param     
    } else if(what=="mergedCounts"){
      query<- lapply(.Object@mergedCounts, function (x) {x[[length(x)]]})     
    } else if(what=="alleleCountsPerBam"){
      query <-  lapply(.Object@alleleCounts, lapply, function (x) {x[[length(x)]]})
    } 
    return(query)
  }
)


#' Report ASB variants
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name BaalChIP.report
#' @description Generates a data.frame per group with significant allele-specific binding (ASB) variants
#' @param .Object An object of the \code{\link{BaalChIP}} class
#' @return A named list, with a data.frame per group. 
#' @details The reported data frame contains the following columns: 
#' \itemize{
#' \item ID: unique identifier string per analysed variant
#' \item Bayes_lower: Lower interval for the estimated allelic ratio (allelic ratio is given by REF / TOTAL)
#' \item Bayes_upper: Uper interval for the estimated allelic ratio (allelic ratio is given by REF / TOTAL)
#' \item Bayes_sig_A: classification of variants into allele-specific towards the reference allele. 1/0 means TRUE/FALSE 
#' \item Bayes_sig_B: classification of variants into allele-specific towards the non-reference allele. 1/0 means TRUE/FALSE 
#' }
#' @seealso \code{\link{summaryASB}}, \code{\link{getASB}}
#' @examples
#' data(BaalObject)
#' report <- BaalChIP.report(BaalObject)
#' 
#' #the reported list is grouped by group_name:
#' names(report)
#' 
#' #check out the report for one of the groups:
#' head(report[["MCF7"]])
#' @export 

setMethod(
  "BaalChIP.report",
  "BaalChIP",
  function(.Object) {
    query<-.Object@ASB 
    return(query)
  }
)


#' Summary of QC 
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name summaryQC
#' @description Generates summary of QC result.
#' @param .Object An object of the \code{\link{BaalChIP}} class
#' @return A list with two elements: 
#' \itemize{
#' \item \code{filtering_stats} containning the number of variants that were filtered out in each filter category and the total number that 'pass' all filters 
#' \item \code{average_stats} containning the average number and average percentage of variants in each filter category, averaged across all analysed groups
#' }
#' @seealso \code{\link{BaalChIP.QC}}, \code{\link{plotQC}}
#' @examples
#'data(BaalObject) 
#'summaryQC(BaalObject)
#' @export 
setMethod(
  "summaryQC",
  "BaalChIP",
  function(.Object) {
    query <- summary_QC(.Object)
    return(query)
  }
)


#' Summary of ASB test
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name summaryASB
#' @description Generates summary of ASB test result.
#' @param .Object An object of the \code{\link{BaalChIP}} class
#' @return A matrix containning the total number of allele-specific variants (TOTAL) and the number of variants allele-specific for the reference (REF) and alternate alleles (ALT).
#' @seealso \code{\link{getASB}}, \code{\link{BaalChIP.report}}
#' @examples
#'data(BaalObject) 
#'summaryASB(BaalObject)
#' @export 
setMethod(
  "summaryASB",
  "BaalChIP",
  function(.Object) {
    query <- summary_ASB(.Object)
    return(query)
  }
)

#' Plots QC results
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name plotQC
#' @description Produces different plots of QC results.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param what A single character value indicating the type of plot. Options: 
#' \itemize{
#' \item \code{overall_pie}: plots the average percentage of variants in each filter category (averaged across all groups analysed)
#' \item \code{boxplot_per_filter}: plots the number of variants that were filtered out per filter category
#' \item \code{barplot_per_group}: plots the number of variants that were filtered out per group
#' }
#' @param addlegend A logical value indicating if legend should be included in the plot (default TRUE)
#' @seealso \code{\link{BaalChIP.QC}}, \code{\link{summaryQC}}
#' #' @examples
#'data(BaalObject) 
#'plotQC(BaalObject, "overall_pie")
#'plotQC(BaalObject, "boxplot_per_filter", addlegend=FALSE)
#'plotQC(BaalObject, "barplot_per_group")
#' @export 
setMethod(
  "plotQC",
  "BaalChIP",
  function(.Object, what= c("overall_pie", "boxplot_per_filter", "barplot_per_group"), addlegend=TRUE) {
    ##-----check input arguments
    if (what == "pie") {what <- "overall_pie"}
    if (what == "boxplot") {what <- "boxplot_per_filter"}
    if (what == "barplot") {what <- "barplot_per_group"}
    BaalChIP.checks(name="plot.what",para=what)
    
    stats <- summary_QC(.Object)
    
    plotfilters(stats=stats, what=what, addlegend=addlegend)
  }
)









