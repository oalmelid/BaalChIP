#BaalChIP: all generics for setMethod
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

##alleleCounts
setGeneric(name="alleleCounts",
                       def=function(.Object,...)
                       {
                               standardGeneric("alleleCounts")
                       }
                       )

##QCfilter
setGeneric(name="QCfilter",
                       def=function(.Object,...)
                       {
                               standardGeneric("QCfilter")
                       }
                       )

##mergeReps
setGeneric(name="filterIntbias",
                       def=function(.Object,...)
                       {
                               standardGeneric("filterIntbias")
                       }
                       )


##mergeReps
setGeneric(name="mergePerGroup",
                       def=function(.Object,...)
                       {
                               standardGeneric("mergePerGroup")
                       }
                       )

##filter1allele
setGeneric(name="filter1allele",
                       def=function(.Object,...)
                       {
                               standardGeneric("filter1allele")
                       }
                       )
                       
##QC step pipeline
setGeneric(name="BaalChIP.QC",
                       def=function(.Object,...)
                       {
                               standardGeneric("BaalChIP.QC")
                       }
                       )

##Detection of ASB variants
setGeneric(name="getASB",
                       def=function(.Object,...)
                       {
                               standardGeneric("getASB")
                       }
                       )

##get slots from BaalChIP
setGeneric(name="BaalChIP.get",
                       def=function(.Object,...)
                       {
                               standardGeneric("BaalChIP.get")
                       }
                       )

##report ASB variants
setGeneric(name="BaalChIP.report",
                       def=function(.Object,...)
                       {
                               standardGeneric("BaalChIP.report")
                       }
                       )

##summaryQC
setGeneric(name="summaryQC",
                       def=function(.Object,...)
                       {
                               standardGeneric("summaryQC")
                       }
                       )

##summaryASB
setGeneric(name="summaryASB",
                       def=function(.Object,...)
                       {
                               standardGeneric("summaryASB")
                       }
                       )


#plotQC
setGeneric(name="plotQC",
                       def=function(.Object,...)
                       {
                               standardGeneric("plotQC")
                       }
                       )








