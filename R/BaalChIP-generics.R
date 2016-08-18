#BaalChIP: all generics for setMethod
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


#' Method alleleCounts
#' @name alleleCounts
#' @rdname alleleCounts
#' @exportMethod alleleCounts
setGeneric(name="alleleCounts",
                       def=function(.Object,
                       min_base_quality=10,
                       min_mapq=15)
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
                       RegionsToFilter = NULL, RegionsToKeep = NULL)
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
                 simulation_script = "local", skipScriptRun=FALSE)
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
                       def=function(.Object, cores=4)
                       {
                               standardGeneric("BaalChIP.run")
                       }
                       )



#' Method getASB
#' @name getASB
#' @rdname getASB
#' @exportMethod getASB
setGeneric(name="getASB",
                       def=function(.Object, Iter = 5000, conf_level = 0.95, cores = 4,
                 RMcorrection = TRUE, RAFcorrection = TRUE)
                       {
                               standardGeneric("getASB")
                       }
                       )


#' Method BaalChIP.get
#' @name BaalChIP.get
#' @rdname BaalChIP.get
#' @exportMethod BaalChIP.get
setGeneric(name="BaalChIP.get",
                       def=function(.Object, what=c("samples","param","alleleCountsPerBam","mergedCounts"))
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




#' BaalObject example dataset. This is a BaalChIP-class object
#' @name BaalObject
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords BaalObject
NULL

#' FAIREexample example dataset. This is a BaalChIP-class object
#' @name FAIREexample
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords FAIREexample
NULL

#' ENCODEexample example dataset. This is a BaalChIP-class object
#' @name ENCODEexample
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords ENCODEexample
NULL


#' Genomic regions of unique mappability
#' @description Unique regions with genomic mappability score equal to 1 selected from DUKE uniqueness mappability track of the UCSC genome browser generated using a window size of 50bp (hg19,  wgEncodeCrgMapabilityAlign50mer table). Used as 'RegionsToKeep' within the QCfilter function so that variants NOT overlapping these regions will be removed.
#' @name UniqueMappability50bp_hg19
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords UniqueMappability50bp_hg19
NULL

#' Genomic regions of collapsed repeats
#' @description Collapsed repeat regions at the 0.1\% threshold (hg19 reference). Used as 'RegionsToFilter' within the QCfilter function so that variants overlapping these regions will be removed.
#' @name pickrell2011cov1_hg19
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @references Pickrell et al., 2011 (\url{http://www.ncbi.nlm.nih.gov/pubmed/21690102})
#' @keywords pickrell2011cov1_hg19
NULL

#' Blacklisted genomic regions
#' @description Blacklisted regions downloaded from the UCSC Genome Browser (mappability track; hg19, wgEncodeDacMapabilityConsensusExcludable and wgEncodeDukeMapabilityRegionsExcludable tables). Used as 'RegionsToFilter' within the QCfilter function so that variants overlapping these regions will be removed.
#' @name blacklist_hg19
#' @docType data
#' @author Ines de Santiago \email{ines.desantiago@@cruk.cam.ac.uk}
#' @keywords blacklist_hg19
NULL
