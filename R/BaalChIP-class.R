##################################################
## BaalChIP.R - Allele specific analysis of  	##
## ChIP-seq data from cancer cell lines         ##
## 23 April 2015                           		  ##
## Ines de Santiago and Wei Liu            		  ##
## ines.desantiago@cruk.cam.ac.uk           	  ##
## Cancer Research UK - University of Cambridge ##
##################################################


#' BaalChIP-class
#' @author Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz
#' @name BaalChIP
#' @param samplesheet A character string indicating the filename for a \code{.csv} file. Column names in the \code{.csv} file should include:
#' \itemize{
#' \item \code{group_name}: identifier string to group samples together
#' \item \code{target}: identifier string for factor (transcription factor, protein)
#' \item \code{replicate_number}: replicate number of sample
#' \item \code{bam_name}: file path for BAM file containing aligned reads for ChIP sample. If duplicated reads are flaged they will not be included in the allelic count data 
#' \item \code{bed_name}: path for BED file containing peaks for ChIP sample 
#' \item \code{SampleID}: identifier string for sample. If not given will use <group_name>_<target>_<replicate_number>
#' }
#' @param hets A named vector with filenames for the \code{.txt} variant files to be used. The names in the vector should correspond to group_name strings in the \code{.csv} samplesheet. Columns names in the \code{.txt} file should include:
#' \itemize{
#' \item ID: unique identifier string per variant. Identifiers have to be unique, and no more than one identifier should be present per data record. If there is no identifier available, then use an arbritary name to name each variant
#' \item CHROM: chromosome identifier from the reference genome per variant (same genome build as BAM and BED files provided)	
#' \item POS: the reference position (1-based)
#' \item REF: reference base. Each base must be one of A,C,G,T in uppercase. Multiple bases are not permitted
#' \item ALT: alternate non-reference base. Each base must be one of A,C,G,T in uppercase. Multiple bases are not permitted
#' \item RAF: [Optional] a value ranging from 0 to 1 for each variant denoting the relative allele frequency (RAF). A value between 0.5 and 1 denotes a bias to the reference allele, and a value between 0 and 0.5 a bias to the alternate allele. If missing, BaalChIP will still run but will not correct for relative allele frequency (copy-number) bias 
#' }
#' @description This S4 class includes a series of methods for detecting allele-specific events from multiple ChIP-seq datasets.
#' @return .Object An object of the \code{\link{BaalChIP}} class.
#' @examples
#'setwd(system.file("test",package="BaalChIP"))
#'samplesheet <- "example.tsv"
#'hets <- c("MCF7"="MCF7_hetSNP.txt", "GM12891"="GM12891_hetSNP.txt")
#'res <- new("BaalChIP", samplesheet=samplesheet, hets=hets)

BaalChIP <- setClass(
	#BaalChIPexperiment class
	"BaalChIP",
	
	slots = c(
                samples = "data.frame",
                hets = "character",
                alleleCounts = "list",
                mergedCounts = "list",
                assayedVar = "list",
                VarTable = "list",
                biasTable = "list",
                ASB = "list",
                param   = "list"
                ),
	
	prototype = c(
		samplesheet = data.frame(),
		hets = "",
		alleleCounts = list(),
		mergedCounts = list(),
		assayedVar = list(),
		VarTable = list(),
		biasTable = list(),
        ASB = list(),
        param = list()
		)
	
)
				




