#BaalChIP: all generics for setMethod
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


#' Method alleleCounts
#' @name alleleCounts
#' @rdname alleleCounts
#' @exportMethod alleleCounts
setGeneric(name="alleleCounts",
                       def=function(.Object,
                       min_base_quality=10,
                       min_mapq=15,
                       verbose=TRUE)
                       {
                               standardGeneric("alleleCounts")
                       }
                       )


#' Method QCfilter
#' @name QCfilter
#' @rdname QCfilter
#' @exportMethod QCfilter
setGeneric(name="QCfilter",
                       def=function(.Object,
                       RegionsToFilter = NULL, RegionsToKeep = NULL,
                       verbose=TRUE)
                       {
                               standardGeneric("QCfilter")
                       }
                       )


#' Method filterIntbias
#' @name filterIntbias
#' @rdname filterIntbias
#' @exportMethod filterIntbias
setGeneric(name="filterIntbias",
                       def=function(.Object, simul_output = NULL, tmpfile_prefix = NULL,
                       simulation_script = "local",
                       alignmentSimulArgs = NULL,
                       skipScriptRun=FALSE,
                       verbose=TRUE)
                       {
                               standardGeneric("filterIntbias")
                       }
                       )



#' Method mergePerGroup
#' @name mergePerGroup
#' @rdname mergePerGroup
#' @exportMethod mergePerGroup
setGeneric(name="mergePerGroup",
                       def=function(.Object)
                       {
                               standardGeneric("mergePerGroup")
                       }
                       )


#' Method filter1allele
#' @name filter1allele
#' @rdname filter1allele
#' @exportMethod filter1allele
setGeneric(name="filter1allele",
                       def=function(.Object)
                       {
                               standardGeneric("filter1allele")
                       }
                       )


#' Method BaalChIP.run
#' @name BaalChIP.run
#' @rdname BaalChIP.run
#' @exportMethod BaalChIP.run
setGeneric(name="BaalChIP.run",
                       def=function(.Object, cores=4, verbose=TRUE)
                       {
                               standardGeneric("BaalChIP.run")
                       }
                       )



#' Method getASB
#' @name getASB
#' @rdname getASB
#' @exportMethod getASB
setGeneric(name="getASB",
                       def=function(.Object, Iter = 5000,
                       conf_level = 0.95, cores = 4,
                       RMcorrection = TRUE,
                       RAFcorrection = TRUE,
                       verbose=TRUE)
                       {
                               standardGeneric("getASB")
                       }
                       )


#' Method BaalChIP.get
#' @name BaalChIP.get
#' @rdname BaalChIP.get
#' @exportMethod BaalChIP.get
setGeneric(name="BaalChIP.get",
                       def=function(.Object, what=c("samples","param","alleleCountsPerBam","mergedCounts","assayedVar", "biasTable"))
                       {
                               standardGeneric("BaalChIP.get")
                       }
                       )


#' Method BaalChIP.report
#' @name BaalChIP.report
#' @rdname BaalChIP.report
#' @exportMethod BaalChIP.report
setGeneric(name="BaalChIP.report",
                       def=function(.Object)
                       {
                               standardGeneric("BaalChIP.report")
                       }
                       )

#' Method summaryQC
#' @name summaryQC
#' @rdname summaryQC
#' @exportMethod summaryQC
setGeneric(name="summaryQC",
                       def=function(.Object)
                       {
                               standardGeneric("summaryQC")
                       }
                       )

#' Method summaryASB
#' @name summaryASB
#' @rdname summaryASB
#' @exportMethod summaryASB
setGeneric(name="summaryASB",
                       def=function(.Object)
                       {
                               standardGeneric("summaryASB")
                       }
                       )


#' Method plotQC
#' @name plotQC
#' @rdname plotQC
#' @exportMethod plotQC
setGeneric(name="plotQC",
                       def=function(.Object, what="barplot_per_group",
                       addlegend=TRUE, plot=TRUE)
                       {
                               standardGeneric("plotQC")
                       }
                       )


#' Method plotSimul
#' @name plotSimul
#' @rdname plotSimul
#' @exportMethod plotSimul
setGeneric(name="plotSimul",
                       def=function(.Object, plot=TRUE)
                       {
                               standardGeneric("plotSimul")
                       }
                       )


#' Method adjustmentBaalPlot
#' @name adjustmentBaalPlot
#' @rdname adjustmentBaalPlot
#' @exportMethod adjustmentBaalPlot
setGeneric(name="adjustmentBaalPlot",
                       def=function(.Object, col = c("green3", "gray50"))
                       {
                               standardGeneric("adjustmentBaalPlot")
                       }
                       )




#' BaalObject example dataset
#' @format A BaalChIP-class object
#' @name BaalObject
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords BaalObject
NULL

#' FAIREexample example dataset
#' @format A BaalChIP-class object
#' @name FAIREexample
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords FAIREexample
#' @source  Downloaded from supplementary data of BaalChIP paper
NULL

#' ENCODEexample example dataset
#' @format A BaalChIP-class object
#' @name ENCODEexample
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @source  Downloaded from supplementary data of BaalChIP paper
#' @keywords ENCODEexample
NULL


#' Genomic regions of unique mappability
#' @description A GRanges object containing unique regions with genomic mappability score equal to 1.
#' Selected from DUKE uniqueness mappability track of the UCSC genome browser (hg19,  wgEncodeCrgMapabilityAlign50mer table).
#'
#' Code used to retrieve these regions:
#' \itemize{
#' \item curl http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw > wgEncodeCrgMapabilityAlign50mer.bw
#' \item ./bigWigToBedGraph wgEncodeCrgMapabilityAlign50mer.bw wgEncodeCrgMapabilityAlign50mer.bedgraph
#' \item awk '{ if ($4 >= 1) print $0 }' wgEncodeCrgMapabilityAlign50mer.bedgraph > wgEncodeCrgMapabilityAlign50mer_UNIQUEregions.bedgraph
#' }
#' Used as 'RegionsToKeep' within the QCfilter function so that variants NOT overlapping these regions will be removed. \cr
#' @details These regions are not applicable to longer reads (> 50bp).
#' @source  Downloaded from \url{http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw}.
#' @name UniqueMappability50bp_hg19
#' @format A GRanges object of 9831690 ranges.
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords UniqueMappability50bp_hg19
NULL

#' Genomic regions of collapsed repeats
#' @description A GRanges object containing collapsed repeat regions at the 0.1\% threshold (hg19 reference).
#' Used as 'RegionsToFilter' within the QCfilter function so that variants overlapping these regions will be removed.
#' @name pickrell2011cov1_hg19
#' @source File available as supplementary data: Pickrell2011_seq.cov1.bed (\url{http://www.ncbi.nlm.nih.gov/pubmed/21690102})
#' @format A GRanges object of 34359 ranges.
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @references Pickrell et al., 2011 (\url{http://www.ncbi.nlm.nih.gov/pubmed/21690102})
#' @keywords pickrell2011cov1_hg19
NULL

#' Blacklisted genomic regions
#' @description A GRanges object containing blacklisted regions identified by the ENCODE and modENCODE consortia.
#' These correspond to artifact regions that tend to show artificially high signal (excessive unstructured anomalous reads mapping).
#' Selected from mappability track of the UCSC genome browser (hg19, wgEncodeDacMapabilityConsensusExcludable and wgEncodeDukeMapabilityRegionsExcludable tables).
#'
#' Code used to retrieve these regions:
#' \itemize{
#' \item curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed > hg19_DACExcludable.txt \cr
#' \item cat hg19_DUKEExcludable.txt hg19_DACExcludable.txt | grep -v "^#" | cut -f 2,3,4,5,6,7 | sort -k1,1 -k2,2n | mergeBed -nms -i stdin > hg19_DUKE_DAC.bed \cr
#' }
#' Used as 'RegionsToFilter' within the QCfilter function so that variants overlapping these regions will be removed.
#' @details Note that these blacklists are applicable to functional genomic data (e.g. ChIP-seq, MNase-seq, DNase-seq, FAIRE-seq) of short reads (20-100bp reads).
#' These are not applicable to RNA-seq or other transcriptome data types.
#' @format A GRanges object of 1378 ranges.
#' @name blacklist_hg19
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords blacklist_hg19
NULL
