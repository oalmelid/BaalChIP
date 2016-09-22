# BaalChIP: all methods


## initialization method
setMethod("initialize", "BaalChIP", function(.Object, samplesheet = NULL, hets = NULL, CorrectWithgDNA = list(),
    .CHECKS = TRUE) {

    ##-----check arguments
    if (missing(samplesheet))
        stop("NOTE: 'samplesheet' is missing!")
    if (missing(hets))
        stop("NOTE: 'hets' is missing!")
    samples <- BaalChIP.checks(name = "samplesheet", samplesheet, .CHECKS = .CHECKS)
    hets <- BaalChIP.checks(name = "hets", list(hets = hets, wd = dirname(samplesheet)), .CHECKS = .CHECKS)
    BaalChIP.checks(name = "gDNA", CorrectWithgDNA, .CHECKS = .CHECKS)

    # passed all checks:
    message("-samplesheet checks: OK!")


    if (is.null(CorrectWithgDNA)) {
        CorrectWithgDNA <- list()
    }

    ##-----check matching names
    cells1 <- unique(samples$group_name)
    checkmatchingnames(names(hets), cells1)
    if (!is.null(names(CorrectWithgDNA))) {
        checkmatchingnames.gDNA(names(CorrectWithgDNA), cells1)
    }

    ##-----initialization
    .Object@samples <- samples
    .Object@hets <- hets
    .Object@gDNA <- CorrectWithgDNA
    .Object@param <- list()
    .Object@simulation_stats <- data.frame()

    .Object
})

setMethod("show", "BaalChIP", function(object) {
    cat(" Type :", class(object), "\n")
    QCstats <- summaryQC(object)[["filtering_stats"]]
    asb <- getBaalSlot(object, "ASB")
    samples <- getBaalSlot(object, "samples")
    cat(" Samples                 :  ", nrow(samples), "\n")
    cat(" Experiments             :  ", unique(samples[, "group_name"]), "\n")
    cat(" Filtering and QC        :  ", ifelse(is.null(QCstats), "None", paste(ncol(QCstats) - 1, "filter(s) applied")),
        "\n")
    cat(" Run allele-specific     :  ", ifelse(length(asb) == 0, "None", "Yes: run BaalChIP.report(object)"),
        "\n")
    cat("\n")
})


#' @export BaalChIP
BaalChIP <- function(samplesheet = NULL, hets = NULL, CorrectWithgDNA = list()) {

    ##-----check arguments
    if (missing(samplesheet))
        stop("NOTE: 'samplesheet' is missing!")
    if (missing(hets))
        stop("NOTE: 'hets' is missing!")
    samples <- BaalChIP.checks(name = "samplesheet", samplesheet, .CHECKS = TRUE)
    hets <- BaalChIP.checks(name = "hets", list(hets = hets, wd = dirname(samplesheet)), .CHECKS = TRUE)
    BaalChIP.checks(name = "gDNA", CorrectWithgDNA, .CHECKS = TRUE)

    if (is.null(CorrectWithgDNA)) {
        CorrectWithgDNA <- list()
    }

    ##-----check matching names
    cells1 <- unique(samples$group_name)
    checkmatchingnames(names(hets), cells1)
    if (!is.null(names(CorrectWithgDNA))) {
        checkmatchingnames.gDNA(names(CorrectWithgDNA), cells1)
    }

    ##-----initialization
    .Object <- new("BaalChIP", samplesheet, hets, CorrectWithgDNA, .CHECKS = TRUE)
    .Object

}

#' Generates allele-specific read count data
#' @import methods
#' @import GenomicRanges
#' @import Rsamtools
#' @import IRanges
#' @importFrom utils read.delim
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom GenomeInfoDb seqlengths
#' @rdname alleleCounts
#' @aliases alleleCounts,BaalChIP-method
#' @author Ines de Santiago
#' @description Generates allele-specific read count data from each BAM ChIP-seq dataset for each variant.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param min_base_quality A numeric value indicating the minimum read base quality below which the base is ignored when summarizing pileup information (default 10).
#' @param min_mapq A numeric value indicating the minimum mapping quality (MAPQ) below which the entire read is ignored (default 15).
#' @param verbose logical. If TRUE reports extra information on the process
#' @note BaalChIP computes allelic counts at each variant position with Rsamtools pileup function. The algorithm follows pileup::Rsamtools by automatically excluding reads flagged as unmapped, secondary, duplicate, or not passing quality controls.
#' @details Utilizes the information within the \code{samples} slot of a BaalChIP object. Will primarily find all variants overlapping peaks. Then, for each variant, computes the number of reads carrying the reference (REF) and alternative (ALT) alleles.
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{alleleCounts} containing a list of GRanges objects.
#' @seealso \code{\link{BaalChIP.get}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'
#'#retrieve alleleCounts:
#'counts <- BaalChIP.get(res, 'alleleCountsPerBam')
#'
#'#alleleCounts are grouped by bam_name and group_name:
#'names(counts)
#'names(counts[['MCF7']])
#'
#'#check out the result for one of the bam files:
#'counts[['MCF7']][[1]]
#' @export
setMethod(f = "alleleCounts", signature = "BaalChIP", function(.Object, min_base_quality = 10, min_mapq = 15,
    verbose = TRUE) {

    ##-----check input arguments
    samples <- getBaalSlot(.Object, "samples")
    hets <- getBaalSlot(.Object, "hets")
    BaalChIP.checks(name = "min_base_quality", min_base_quality)
    BaalChIP.checks(name = "min_mapq", min_mapq)

    ##-----assign parameters
    .Object@param$QCparam <- list(min_base_quality = min_base_quality, min_mapq = min_mapq)

    ##-----compute allele counts
    res_per_bam <- applyAlleleCountsPerBam(samples, hets, min_base_quality, min_mapq, verbose = verbose)
    .Object@alleleCounts <- res_per_bam

    return(.Object)
})


#' Removes variants that may be problematic for identification of allele-specific events
#' @import methods
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @import GenomicRanges
#' @import IRanges
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels
#' @author Ines de Santiago
#' @rdname QCfilter
#' @aliases QCfilter,BaalChIP-method
#' @description Quality control step for removing variants that may be problematic for identification of allele-specific events.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param RegionsToFilter a named list of GRanges objects. Variants overlapping these regions will be removed.
#' @param RegionsToKeep a named list of GRanges objects. Works in an oposite way to 'RegionstoFilter', variants NOT overlapping these regions will be removed.
#' @param verbose logical. If TRUE reports extra information on the process
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{alleleCounts} containing a list of GRanges objects that pass filters.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data('blacklist_hg19')
#'data('pickrell2011cov1_hg19')
#'data('UniqueMappability50bp_hg19')
#'res <- QCfilter(res,
#'                RegionsToFilter=list('blacklist'=blacklist_hg19,
#'                'highcoverage'=pickrell2011cov1_hg19),
#'                RegionsToKeep=list('UniqueMappability'=UniqueMappability50bp_hg19))
#'
#'#check results
#'plotQC(res,'barplot')
#'summaryQC(res)
#' @export
setMethod(f = "QCfilter", signature = "BaalChIP", function(.Object, RegionsToFilter = NULL, RegionsToKeep = NULL,
    verbose = TRUE) {

    ##-----check input arguments
    RegionsToFilter <- BaalChIP.checks(name = "RegionsToFilter", RegionsToFilter)
    RegionsToKeep <- BaalChIP.checks(name = "RegionsToKeep", RegionsToKeep)

    ##-----check that mergedCounts and ASB have not been computed yet
    asb <- getBaalSlot(.Object, "ASB")
    merged <- .Object@mergedCounts
    if (!all(length(asb) == 0 & length(merged) == 0)) {
        # message(prompt='Running this QC step at the single BAM files level will delete any downstream
        # analysis. Continue?[y/n]') n <- scan('stdin', character(), n=1) if (n == 'y') {
        .Object@ASB <- list()
        .Object@mergedCounts <- list()
        .Object@assayedVar <- list()
        .Object@biasTable <- list()
        .Object@VarTable <- list()
        # }else{stop('Interrupted by the user')}
    }

    ##-----assign parameters
    counts_per_bam <- getBaalSlot(.Object, "alleleCounts")
    .Object@param$QCparam$RegionsToFilter <- names(RegionsToFilter)
    .Object@param$QCparam$RegionsToKeep <- names(RegionsToKeep)

    ##-----do not run if alleleCounts not found
    if (length(counts_per_bam) == 0) {
        stop("Please run alleleCounts function before running QCfilter")
    }

    ##-----apply filters
    res_per_bam <- applyFiltersPerBam(counts_per_bam, RegionsToFilter, RegionsToKeep, verbose = verbose)
    .Object@alleleCounts <- res_per_bam

    return(.Object)
})


#' Removes variants for which simulated single-end reads don't align correctly to the original simulated position.
#' @import methods
#' @importFrom utils read.delim
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom utils write.table
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomicAlignments
#' @author Ines de Santiago
#' @rdname filterIntbias
#' @aliases filterIntbias,BaalChIP-method
#' @description Filters the data frame available within a \code{\link{BaalChIP}} object (slot \code{alleleCounts}). This filter performs simulations of reads of the same length as the original ChIP-seq reads, aligns the simulated reads to the genome, calculates the allelic ratios for each variant and finally ignores those variants for which the allelic ratio (REF/TOTAL) is different than 0.5.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param simul_output a non-empty character vector giving the directory of where to save the FASTQ and BAM files generated by the function. If NULL, a random directory under the current working directory will be generated.
#' @param tmpfile_prefix an optional character vector giving the initial part of the name of the FASTQ and BAM files generated by the function. If NULL, a random name will be generated.
#' @param simulation_script the file path for simulation script containing the instructions of simulation and alignment commands. If NULL, the default simulation script distributed with BaalChIP ('extra/simulation_run.sh') will be used.
#' @param skipScriptRun a logical value indicating if simulation BAM files should not be generated. If TRUE BaalChIP will look for the BAM files in the 'simul_output/temp_prefix' (default is FALSE).
#' @param alignmentSimulArgs a character vector with arguments to the simulation script. If NULL no arguments are passed.
#' @param verbose logical. If TRUE reports extra information on the process
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{alleleCounts} containing a list of GRanges objects that pass filters.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'skipScriptRun=TRUE #For demonstration purposes only (read details in vignette)
#'res <- filterIntbias(res,
#'       simul_output=system.file('test/simuloutput',package='BaalChIP'),
#'       tmpfile_prefix='c67c6ec6c433', skipScriptRun=TRUE)
#'
#'#check results
#'plotSimul(res)
#'summaryQC(res)
#' @export
setMethod(f = "filterIntbias", signature = "BaalChIP", function(.Object, simul_output = NULL, tmpfile_prefix = NULL,
    simulation_script = "local", alignmentSimulArgs = NULL, skipScriptRun = FALSE, verbose = TRUE) {

    ##-----check input arguments
    simul_output <- BaalChIP.checks(name = "simul_output", list(dir = simul_output, prefix = tmpfile_prefix))
    BaalChIP.checks(simulation_script, "simulation_script")
    BaalChIP.checks(skipScriptRun, "skipScriptRun")

    ##-----assign parameters
    res_per_bam <- getBaalSlot(.Object, "alleleCounts")
    samples <- getBaalSlot(.Object, "samples")
    .Object@param$QCparam$FiltIntBias <- TRUE
    .Object@param$QCparam$simul_output <- simul_output

    ##-----do not run if alleleCounts not found
    if (length(res_per_bam) == 0) {
        stop("Please run alleleCounts (and optionally QCfilter) functions before running filterIntbias")
    }

    ##-----read length per sample
    samples$readlen <- applyReadlenPerBam(samples, res_per_bam, verbose = verbose)

    ##---- run simulations
    # DEBUG if (testingScript) { #save so far before running simulations tmpname <-
    # file.path(dirname(samplesheet),paste0(tmpfile_prefix,'_step1.Rda')) save(samples, res_per_bam,
    # file = tmpname) }

    # simulate existsTest <- exists('testingScript') #if does not exist will return FALSE if (!
    # existsTest) {testingScript = FALSE} if (!is.logical(testingScript)) {stop('cannot continue,
    # testingScript has to be either TRUE or FALSE')}
    testingScript <- skipScriptRun
    simures <- applySim(samples, res_per_bam, simul_output = simul_output, simulation_script = simulation_script,
        alignmentSimulArgs = alignmentSimulArgs, testingScript = testingScript, verbose = verbose)

    ##----apply filter
    res_per_bam <- applyIntBiasFilterPerBam(samples, res_per_bam, simcounts = simures[["simcounts"]],
        verbose = verbose)

    # save so far DEBUG if (testingScript) { tmpname <-
    # file.path(dirname(samplesheet),paste0(tmpfile_prefix,'_step2.Rda')) save(samples, res_per_bam,
    # simures, file = tmpname) }

    .Object@simulation_stats <- simures[["simulation_stats"]]
    .Object@alleleCounts <- res_per_bam
    return(.Object)
})


#' Merges allele-specific read count data per group
#' @import methods
#' @import GenomicRanges
#' @author Ines de Santiago
#' @rdname mergePerGroup
#' @aliases mergePerGroup,BaalChIP-method
#' @description Merges all ChIP-seq datasets within a group of samples creating a data.frame that contains allele-specific read count data for all variants that need to be analysed.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @details if QCfilter has been applied, will use the most up-to-date variant set available for each individual BAM file (after QC). Missing values are allowed for heterozygous variants that are not available (e.g. do not pass filter for a particular ChIP-seq dataset).
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{mergedCounts} containing a data.frame of merged samples per group.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data('blacklist_hg19')
#'data('pickrell2011cov1_hg19')
#'data('UniqueMappability50bp_hg19')
#'res <- QCfilter(res,
#'                RegionsToFilter=list('blacklist'=blacklist_hg19,
#'                'highcoverage'=pickrell2011cov1_hg19),
#'                RegionsToKeep=list('UniqueMappability'=UniqueMappability50bp_hg19))
#'
#'res <- mergePerGroup(res)
#'
#'#retrieve mergedCounts:
#'counts <- BaalChIP.get(res, 'mergedCounts')
#'
#'#mergedCounts are grouped by group_name:
#'names(counts)
#'sapply(counts, dim)
#'
#'#check out the result for one of the groups:
#'head(counts[[1]])
#' @export
setMethod(f = "mergePerGroup", signature = "BaalChIP", function(.Object) {

    ##-----assign parameters
    res_per_bam <- getBaalSlot(.Object, "alleleCounts")
    samples <- getBaalSlot(.Object, "samples")

    ##-----do not run if alleleCounts not found
    if (length(res_per_bam) == 0) {
        stop("Please run alleleCounts (and optionally QCfilter and/or filtIntBias) functions before running mergePerGroup")
    }

    ##-----apply filters
    res_merged <- applyMergeResults(samples, res_per_bam)
    .Object@mergedCounts <- res_merged

    return(.Object)
})

#' Filters out variants with only 1 observed allele
#' @import methods
#' @import GenomicRanges
#' @author Ines de Santiago
#' @rdname filter1allele
#' @aliases filter1allele,BaalChIP-method
#' @description Filters the data frame available within a \code{\link{BaalChIP}} object (slot \code{mergedCounts}). This filter ignores variants for which only one allele is observed after pooling ChIP-seq reads from all datasets.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{mergedCounts} containing a data.frame of merged samples per group with variants that pass the filter.
#' @seealso \code{\link{BaalChIP.get}}, \code{\link{plotQC}}, \code{\link{summaryQC}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'data('blacklist_hg19')
#'data('pickrell2011cov1_hg19')
#'data('UniqueMappability50bp_hg19')
#'res <- QCfilter(res,
#'                RegionsToFilter=list('blacklist'=blacklist_hg19,
#'                'highcoverage'=pickrell2011cov1_hg19),
#'                RegionsToKeep=list('UniqueMappability'=UniqueMappability50bp_hg19))
#'
#'res <- mergePerGroup(res)
#'res <- filter1allele(res)
#'
#'#retrieve mergedCounts:
#'counts <- BaalChIP.get(res, 'mergedCounts')
#'
#'#mergedCounts are grouped by group_name:
#'names(counts)
#'sapply(counts, dim)
#'
#'#check out the result for one of the groups:
#'head(counts[[1]])
#' @export
setMethod(f = "filter1allele", signature = "BaalChIP", function(.Object) {

    ##-----assign parameters
    res_merged <- .Object@mergedCounts

    ##-----do not run if res_merged not found
    if (length(res_merged) == 0) {
        stop("Please run mergePerGroup function before running filter1allele")
    }

    ##-----apply filters
    res_merged <- applyFilter1allele(res_merged)
    .Object@mergedCounts <- res_merged

    return(.Object)
})

#' BaalChIP pipeline - allele counts,  QC and getting ASB variants
#' @import methods
#' @importFrom utils data
#' @importFrom utils read.delim
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom utils write.table
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomicAlignments
#' @import IRanges
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats dbeta
#' @importFrom coda as.mcmc
#' @import foreach
#' @importFrom coda HPDinterval
#' @author Ines de Santiago
#' @rdname BaalChIP.run
#' @aliases BaalChIP.run,BaalChIP-method
#' @description BaalChIP.run is a wrapper convenience function, to compute allele counts and perform quality controls in one step. This function will use the package's defaults.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param cores number of cores for parallel computing (default is 4).
#' @param verbose logical. If TRUE reports extra information on the process
#' @details This function is a wrapper of the following functions: \code{\link{alleleCounts}}, \code{\link{QCfilter}}, \code{\link{mergePerGroup}}, \code{\link{filter1allele}}, \code{\link{getASB}}
#' @return An object of the \code{\link{BaalChIP}} class.
#' @seealso \code{\link{summaryQC}}, \code{\link{plotQC}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- BaalChIP.run(res, cores=2)
#'
#'#summary of the QC step
#'summaryQC(res)
#'#summary of the ASB step
#'summaryASB(res)
#' @export
setMethod(f = "BaalChIP.run", signature = "BaalChIP", function(.Object, cores = 4, verbose = TRUE) {
    # simul_output=NULL)

    ##-----check input arguments
    BaalChIP.checks(name = "cores", cores)

    ##----- compute allele counts
    .Object <- alleleCounts(.Object, min_base_quality = 10, min_mapq = 15, verbose = verbose)

    ##-----run QC step
    data("blacklist_hg19")
    data("pickrell2011cov1_hg19")
    data("UniqueMappability50bp_hg19")
    .Object <- QCfilter(.Object, RegionsToFilter = c("blacklist_hg19", "pickrell2011cov1_hg19"), RegionsToKeep = c("UniqueMappability50bp_hg19"),
        verbose = verbose)

    # if (FiltIntBias == FALSE & !is.null(simul_output)) { warning (paste('will not use 'simul_output'
    # because FiltIntBias is FALSE'))

    ##-----merge replicates
    if (verbose) {
        message("-merging replicated samples...")
    }
    .Object <- mergePerGroup(.Object)

    ##-----filter 'Only1Allele'
    if (verbose) {
        message("-filtering out SNPs with only 1 observed allele...")
    }
    .Object <- filter1allele(.Object)

    if (verbose) {
        message("-QC complete!\n")
    }

    ##-----get ASB
    if (verbose) {
        message("-geting ASB counts...")
    }
    .Object <- getASB(.Object, Iter = 5000, conf_level = 0.95, cores = cores, verbose = verbose)

    .Object
})

#' Identifies allele-specific binding events
#' @import methods
#' @importFrom utils read.delim
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats dbeta
#' @importFrom coda as.mcmc
#' @import foreach
#' @importFrom coda HPDinterval
#' @author Wei Liu, Ke Yuan, Ines de Santiago
#' @rdname getASB
#' @aliases getASB,BaalChIP-method
#' @description getASB identifies allele-specific binding events using a bayesian framework.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param Iter Maximum number of iterations (default 5000).
#' @param conf_level Confidence interval in the estimated allelic ratio (default 0.95).
#' @param cores number of cores for parallel computing (default is 4).
#' @param RMcorrection Logical value indicating if reference mapping (RM) bias should be applied (default TRUE). If FALSE will not correct for reference allele mapping bias. If TRUE will estimate the RM bias from the overall reference allele proportion.
#' @param RAFcorrection Logical value indicating if relative allele frequency (RAF) bias correction should be applied (default TRUE). If TRUE will read RAF values for each variant from \code{hets} files (RAF column name). If FALSE will not correct for relative allele frequency bias.
#' @param verbose logical. If TRUE reports extra information on the process
#' @return An updated \code{\link{BaalChIP}} object with the slot \code{ASB} containing variants identified as allele-specific.
#' @seealso \code{\link{summaryASB}}, \code{\link{BaalChIP.report}}
#' @examples
#'setwd(system.file('test',package='BaalChIP'))
#'samplesheet <- 'exampleChIP.tsv'
#'hets <- c('MCF7'='MCF7_hetSNP.txt', 'GM12891'='GM12891_hetSNP.txt')
#'res <- BaalChIP(samplesheet=samplesheet, hets=hets)
#'res <- alleleCounts(res, min_base_quality=10, min_mapq=15)
#'res <- mergePerGroup(res)
#'res <- getASB(res, cores=2)
#'
#'#summary - number of significant ASB variants
#'summaryASB(res)
#'
#'#report result
#'res <- BaalChIP.report(res)
#' @export
setMethod("getASB", "BaalChIP", function(.Object, Iter = 5000, conf_level = 0.95, cores = 4, RMcorrection = TRUE,
    RAFcorrection = TRUE, verbose = TRUE) {

    ##-----check input arguments
    BaalChIP.checks(name = "Iter", Iter)
    BaalChIP.checks(name = "conf_level", conf_level)
    BaalChIP.checks(name = "cores", cores)
    BaalChIP.checks(name = "RMcorrection", RMcorrection)
    BaalChIP.checks(name = "RAFcorrection", RAFcorrection)

    ##-----update object with 'assayedVar'
    assayedVar <- BaalChIP.get(.Object, "mergedCounts")  #last filtered mergedCounts
    .Object@assayedVar <- assayedVar

    if (length(assayedVar) == 0) {
        stop("there are no variants in 'mergedCounts' slot. Make sure you have run alleleCounts and mergePerGroup functions before continuing")
    }

    #------get VarTable
    VarTable <- addVarTable(.Object, RAFcorrection = RAFcorrection, verbose = verbose)
    .Object@VarTable <- VarTable

    ##-----check matching cellnames
    samples <- getBaalSlot(.Object, "samples")
    assayedVar <- getBaalSlot(.Object, "assayedVar")
    VarTable <- getBaalSlot(.Object, "VarTable")
    cells1 <- unique(samples$group_name)
    checkmatchingnames(names(assayedVar), cells1)
    checkmatchingnames(names(VarTable), cells1)

    ##-----run
    Expnames <- names(assayedVar)
    results <- list()
    biasTable <- list()
    applyedCorrection <- list()

    for (ID in Expnames) {
        message("... calculating ASB for: ", ID)
        assayed <- assayedVar[[ID]]
        GTtable <- VarTable[[ID]]

        # if no RAF correction replace any RAF values by 0.5 (even if they were given as input args..)
        if (!RAFcorrection) {
            GTtable$RAF <- 0.5
        }

        # get bias table (variants ordered equally between counts table and biastable)
        if (RMcorrection) {
            ARestimate <- estimateRefBias(assayed, GTtable, min_n = 200)
        } else {
            ARestimate = NULL
        }

        # get biasTable
        result <- getbiasTable(assayed, GTtable, ARestimate)
        counts <- result[[1]]
        biastable <- result[[2]]
        biasparam <- getbiasparam(biastable)

        # run bayes
        if (nrow(counts) > 0) {
            Bayes_report <- runBayes(counts = counts, bias = biastable, Iter = Iter, conf_level = conf_level,
                cores = cores)
        } else {
            message("no variants left for ", ID)
            Bayes_report <- data.frame()
        }

        # append results
        results[[ID]] <- Bayes_report
        biasTable[[ID]] <- biastable
        applyedCorrection[[ID]] <- biasparam
    }


    ##-----assign parameters
    applyedCorrection <- t(do.call("rbind", applyedCorrection))
    .Object@param$ASBparam <- list(Iter = Iter, conf_level = conf_level, applyedCorrection = applyedCorrection)

    ##-----update status and return
    if (!all(sapply(results, nrow) == 0)) {
        .Object@biasTable <- biasTable
        .Object@ASB <- results
    }
    message("-ASB identification complete!")
    .Object
})


#' Get slots from a BaalChIP object
#' @import methods
#' @author Ines de Santiago
#' @rdname BaalChIP.get
#' @aliases BaalChIP.get,BaalChIP-method
#' @description Get information from individual slots in a BaalChIP object.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param what a single character value specifying which information should be retrieved. Options: 'samples', 'param', 'alleleCountsPerBam', 'mergedCounts', 'assayedVar', 'biasTable'.
#' @return The slot content from an object of the \code{\link{BaalChIP}} class.
#' @examples
#'data('BaalObject')
#'
#'#samples data spreadsheet and hets:
#'BaalChIP.get(BaalObject,what='samples')
#'
#'#parameters used within run:
#'BaalChIP.get(BaalObject,what='param')
#'
#'#retrieve a GRanges list with allele-specific read counts per BAM file:
#'counts <- BaalChIP.get(BaalObject,what='alleleCountsPerBam')
#'counts[['MCF7']][[1]]
#'
#'#retrieve a data.frame with allele-specific read counts per group:
#'counts <- BaalChIP.get(BaalObject,what='mergedCounts')
#'head(counts[[1]])
#' @export
setMethod("BaalChIP.get", "BaalChIP", function(.Object, what = c("samples", "param", "alleleCountsPerBam",
    "mergedCounts", "assayedVar")) {
    ##-----check input arguments
    BaalChIP.checks(name = "get.what", param = what)
    ##-----get query
    query <- NULL
    if (what == "samples") {
        query <- list(samples = .Object@samples, hets = .Object@hets)
    } else if (what == "param") {
        query <- .Object@param
    } else if (what == "mergedCounts") {
        query <- lapply(.Object@mergedCounts, function(x) {
            x[[length(x)]]
        })
    } else if (what == "alleleCountsPerBam") {
        query <- lapply(.Object@alleleCounts, lapply, function(x) {
            x[[length(x)]]
        })
    } else if (what == "assayedVar") {
        query <- .Object@assayedVar
    } else if (what == "biasTable") {
        query <- .Object@biasTable
    }
    return(query)
})

getBaalSlot <- function(.Object, what) {
    match.arg(what, c("samples", "param", "mergedCounts", "alleleCounts", "assayedVar", "biasTable",
        "hets", "simulation_stats", "ASB", "VarTable", "gDNA"))
    switch(what, samples = {
        query <- .Object@samples
    }, param = {
        query <- .Object@param
    }, mergedCounts = {
        query <- .Object@mergedCounts
    }, alleleCounts = {
        query <- .Object@alleleCounts
    }, assayedVar = {
        query <- .Object@assayedVar
    }, biasTable = {
        query <- .Object@biasTable
    }, hets = {
        query <- .Object@hets
    }, simulation_stats = {
        query <- .Object@simulation_stats
    }, ASB = {
        query <- .Object@ASB
    }, VarTable = {
        query <- .Object@VarTable
    }, gDNA = {
        query <- .Object@gDNA
    }, )
    return(query)
}

#' Report ASB variants
#' @import methods
#' @author Ines de Santiago
#' @rdname BaalChIP.report
#' @aliases BaalChIP.report,BaalChIP-method
#' @description Generates a data.frame per group with all variants and a label for all identified allele-specific binding (ASB) variants.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @return A named list, with a data.frame per group.
#' @details The reported data frame contains the following columns:
#' \itemize{
#' \item ID: unique identifier string per analysed variant.
#' \item CHROM: chromosome identifier from the reference genome per variant.
#' \item POS: the reference position (1-based).
#' \item REF: reference base. Each base must be one of A,C,G,T in uppercase.
#' \item ALT: alternate non-reference base. Each base must be one of A,C,G,T in uppercase.
#' \item REF.counts: pooled counts of all reads with the reference allele.
#' \item ALT.counts: pooled counts of all reads with the non-reference allele.
#' \item Total.counts: pooled counts of all reads (REF + ALT).
#' \item AR: allelic ratio calculated directly from sequencing reads (REF / TOTAL).
#' \item RMbias: numerical value indicating the value estimated and applied by BaalChIP for the reference mapping bias. A value between 0.5 and 1 denotes a bias to the reference allele, and a value between 0 and 0.5 a bias to the alternative allele.
#' \item RAF: numerical value indicating the value applied by BaalChIP for the relative allele frequency (RAF) bias correction. A value between 0.5 and 1 denotes a bias to the reference allele, and a value between 0 and 0.5 a bias to the alternative allele.
#' \item Bayes_lower: lower interval for the estimated allelic ratio (allelic ratio is given by REF / TOTAL).
#' \item Bayes_upper: upper interval for the estimated allelic ratio (allelic ratio is given by REF / TOTAL).
#' \item Corrected.AR: average estimated allelic ratio (average between Bayes_lower and Bayes_upper). A value between 0.5 and 1 denotes a bias to the reference allele, and a value between 0 and 0.5 a bias to the alternative allele.
#' \item isASB: logical value indicating BaalChIP's classification of variants into allele-specific.
#' }
#' @seealso \code{\link{summaryASB}}, \code{\link{getASB}}
#' @examples
#'data('BaalObject')
#'report <- BaalChIP.report(BaalObject)
#'
#'#the reported list is grouped by group_name:
#'names(report)
#'
#'#check out the report for one of the groups:
#'head(report[['MCF7']])
#' @export
setMethod("BaalChIP.report", "BaalChIP", function(.Object) {
    samples <- getBaalSlot(.Object, "samples")
    group_names <- unique(samples[["group_name"]])
    query <- lapply(group_names, function(x) getReport(.Object, group_name = x))
    names(query) <- group_names
    if (all(sapply(query, is.null))) {
        query <- NULL
    }
    return(query)
})

#' Summary of QC
#' @import methods
#' @importFrom reshape2 melt
#' @importFrom doBy summaryBy
#' @author Ines de Santiago
#' @rdname summaryQC
#' @aliases summaryQC,BaalChIP-method
#' @description Generates summary of QC result.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @return A list with two elements:
#' \itemize{
#' \item \code{filtering_stats} containing the number of variants that were filtered out in each filter category and the total number that 'pass' all filters.
#' \item \code{average_stats} containing the average number and average percentage of variants in each filter category, averaged across all analysed groups.
#' }
#' @seealso \code{\link{BaalChIP.run}}, \code{\link{plotQC}}
#' @examples
#'data('BaalObject')
#'summaryQC(BaalObject)
#' @export
setMethod("summaryQC", "BaalChIP", function(.Object) {
    query <- summary_QC(.Object)
    if (!is.null(query)) {
        return(query)
    } else {
        return(invisible(NULL))
    }
})


#' Summary of ASB test
#' @import methods
#' @author Ines de Santiago
#' @rdname summaryASB
#' @aliases summaryASB,BaalChIP-method
#' @description Generates summary of ASB test result.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @return A matrix containing the total number of allele-specific variants (TOTAL) and the number of variants allele-specific for the reference (REF) and alternate alleles (ALT).
#' @seealso \code{\link{getASB}}, \code{\link{BaalChIP.report}}
#' @examples
#'data('BaalObject')
#'summaryASB(BaalObject)
#' @export
setMethod("summaryASB", "BaalChIP", function(.Object) {
    query <- summary_ASB(.Object)
    return(query)
})

#' Plots QC results
#' @import methods
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom doBy summaryBy
#' @importFrom scales alpha
#' @importFrom graphics plot
#' @author Ines de Santiago
#' @rdname plotQC
#' @aliases plotQC,BaalChIP-method
#' @description Produces different plots of QC results.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param what A single character value indicating the type of plot. Options:
#' \itemize{
#' \item \code{overall_pie}: plots the average percentage of variants in each filter category (averaged across all groups analysed).
#' \item \code{boxplot_per_filter}: plots the number of variants that were filtered out per filter category.
#' \item \code{barplot_per_group}: plots the number of variants that were filtered out per group.
#' }
#' @param addlegend A logical value indicating if legend should be included in the plot (default TRUE).
#' @param plot a logical value to whether it should plot (TRUE) or not (FALSE). Default is TRUE.
#' @return A plot
#' @seealso \code{\link{BaalChIP.run}}, \code{\link{summaryQC}}
#' @examples
#'data('BaalObject')
#'plotQC(BaalObject, what = 'overall_pie')
#'plotQC(BaalObject, what =  'boxplot_per_filter', addlegend=FALSE)
#'plotQC(BaalObject, what =  'barplot_per_group')
#' @export
setMethod("plotQC", "BaalChIP", function(.Object, what = "barplot_per_group", addlegend = TRUE, plot = TRUE) {
    ##-----check input arguments
    match.arg(what, c("overall_pie", "boxplot_per_filter", "barplot_per_group", "pie", "boxplot", "barplot"))
    switch(what, pie = {
        what <- "overall_pie"
    }, boxplot = {
        what <- "boxplot_per_filter"
    }, barplot = {
        what <- "barplot_per_group"
    })
    BaalChIP.checks(name = "plot.what", param = what)
    stats <- summary_QC(.Object)
    p <- plotfilters(stats = stats, what = what, addlegend = addlegend, plot = plot)
    if (!plot) {
        return(p)
    }
})

#' Plots the percentage of correctly aligned simulated reads
#' @import methods
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom doBy summaryBy
#' @importFrom scales alpha
#' @importFrom graphics plot
#' @author Ines de Santiago
#' @rdname plotSimul
#' @aliases plotSimul,BaalChIP-method
#' @description Produces a plot of the proportion of variants that displayed the correct number of mapped simulated reads.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param plot a logical value to whether it should plot (TRUE) or not (FALSE). Default is TRUE.
#' @return A plot
#' @examples
#'data('BaalObject')
#'plotSimul(BaalObject)
#' @export
setMethod("plotSimul", "BaalChIP", function(.Object, plot = TRUE) {

    simulation_stats <- getBaalSlot(.Object, "simulation_stats")
    p <- plot.simul(simulation_stats, plot = plot)
    if (!plot) {
        return(p)
    }
})

#' Plots Allelic Ratios before and after adjustment
#' @import methods
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom doBy summaryBy
#' @importFrom scales alpha
#' @importFrom graphics plot
#' @author Ines de Santiago
#' @rdname adjustmentBaalPlot
#' @aliases adjustmentBaalPlot,BaalChIP-method
#' @description Produces a density plot of the distribution of allelic ratios (REF/TOTAL) before and after BaalChIP adjustment for RM and RAF biases.
#' @param .Object An object of the \code{\link{BaalChIP}} class.
#' @param col A character vector indicating the colours for the density plot ( default is c( 'green3','gray50') ).
#' @return A plot
#' @seealso \code{\link{BaalChIP.report}}, \code{\link{summaryQC}}
#' @examples
#'data('BaalObject')
#'adjustmentBaalPlot(BaalObject)
#'adjustmentBaalPlot(BaalObject, col=c('blue','pink'))
#' @export
setMethod("adjustmentBaalPlot", "BaalChIP", function(.Object, col = c("green3", "gray50")) {

    query <- BaalChIP.report(.Object)
    if (!is.null(query)) {
        plotadjustment(query, col = col)
    } else {
        message("nothing to plot")
    }
})







