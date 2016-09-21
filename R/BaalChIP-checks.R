#BaalChIP: functions for argument checking
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

checkFileExists <- function(fname, wd) {
    #checks if fname exists, if not will try to see inside wd. Returns right fname
    if (file.exists(fname)) {
        return(fname)
    }else{
        fname2 <- file.path(wd, fname)
        if (file.exists(fname2)){
            return(fname2)
        }else
            stop(paste('file does not exist:', fname),call.=FALSE)
        }
}

readsamplesheet <- function(samplesheet, .CHECKS=TRUE) {
    # Read samplesheet
    samples <- read.delim(samplesheet,stringsAsFactors=FALSE)

    #check if colnames exist
    colnames_in_file = c('group_name', 'target', 'replicate_number', 'bam_name', 'bed_name')
    if (! all(colnames_in_file %in% colnames(samples))) {
            stop(paste('samplesheet must contain the following colnames:',  paste(colnames_in_file, collapse=",")),call.=FALSE)
    }

    # make sampleIDs
    if (! "sampleID" %in% colnames(samples)){
      sampleID <- paste(samples[,'group_name'],samples[,'target'],samples[,'replicate_number'], sep='_')
        samples$sampleID <- sampleID
    }

    #Check duplicated names
    if (any(duplicated(sampleID))) {
      i <- which(duplicated(sampleID))
      stop(paste('Error: duplicated sample IDs', sampleID[i]),call.=FALSE)
    }

    #check if all files exist in samplesheet
    if (.CHECKS) {
        wd <- dirname(samplesheet)
        samples$bam_name <- sapply(samples$bam_name, checkFileExists, wd = wd)
        a <- sapply(paste0(samples$bam_name,".bai"), checkFileExists, wd = dirname(samples$bam_name))
        samples$bed_name <- sapply(samples$bed_name, checkFileExists, wd = wd)
    }

    return(samples)
}

readhettables <- function(hets, wd, .CHECKS=TRUE) {
    if (.CHECKS) {
        hets <- sapply(hets, checkFileExists, wd = wd)
    }
    return(hets)
}

readbams <- function(bamlist, .CHECKS=TRUE) {
    if (.CHECKS & length(bamlist)>0) {
        for (filename in unlist(bamlist)) {
            if (!file.exists(filename)) { stop(paste('bam file does not exist:',filename),call.=FALSE) }
        }
    }
}

trycreatedir <- function(dir) {
    if (!is.null(dir)) {
        if (!file.exists(dir)){
            try(dir.create(dir), silent=TRUE)
            if (!file.exists(dir)){
                stop (paste('failed to create simul_output dir:', dir),call.=FALSE)
            }
        }
    }
}

createtempdir <- function(dirlist) {

    dir1 <- dirlist[["dir"]]
    prefix <- dirlist[["prefix"]]

    #create output directory
    if (is.null(dir1)) {dir1=tempfile(tmpdir=getwd(), pattern="BaalChIP_simulOUT_")}
    trycreatedir(dir1)

    #add the prefix
    if (!is.null(prefix)) { if(prefix == "") {prefix <- NULL}}
    if (is.null(prefix)) {
        dir1 <- tempfile(tmpdir=dir1, pattern="")
    }else{
        dir1 <- file.path(dir1, prefix)
    }
    return(dir1)
}


checkmatchingnames <- function(names1, names2) {

    if (is.null(names2)) {stop("cannot check names")}

    #will give a warning for unmatching names
    if (!all(names1 %in% names2)) {warning("group name not found in samplesheet: ", paste(setdiff(names1, names2), collapse=","), call.=FALSE)}

    #will Fail if all names2 not in names1
    if (!all(names2 %in% names1)) {stop("samplesheet group name not found in 'hets' vector: ", paste(setdiff(names2, names1), collapse=","), call.=FALSE)}

}

checkmatchingnames.gDNA <- function(gDNA_names, samplesheet_names) {
    if (!all(gDNA_names %in% samplesheet_names)) {warning("gDNA group name not found in samplesheet: ", paste(setdiff(gDNA_names, samplesheet_names), collapse=","), call.=FALSE)}

}

BaalChIP.checks <- function(name, param, .CHECKS= TRUE){

    if (name == "samplesheet")
        {
        samples <- readsamplesheet(param, .CHECKS = .CHECKS)
        return(samples)
    }

    if (name == "hets")
        {
        hets <- readhettables(param[["hets"]], param[["wd"]], .CHECKS = .CHECKS)
        return(hets)
    }

    if (name == "gDNA")
        {
        readbams(param, .CHECKS = .CHECKS)
    }

    if (name == "min_base_quality") {
        if (class(param) != 'numeric') {
            stop ('min_base_quality must be a numberic value',call.=FALSE)
        }
    }

    if (name == "min_mapq") {
        if (class(param) != 'numeric') {
            stop ('min_mapq must be a numberic value',call.=FALSE)
        }
    }

    if (name == "RegionsToFilter" | name == "RegionsToKeep") {
        if (is.character(param)) {
            region_names <- param #Regions contains the names to lists in the environment
            param <- lapply(region_names, get)
            names(param) <- region_names
        }

        if (!is.null(param) &  !is.character(param) & is.null(names(param)))
            {
            stop("Error: RegionsToFilter needs to be a named list. No names found",call.=FALSE)
        }

        return(param)
    }

    if (name == "FiltIntBias") {
        if (!is.logical(param)) {
            stop (paste('FiltIntBias must be either TRUE or FALSE'),call.=FALSE)
        }
    }

    if (name == "Only1Allele") {
        if (!is.logical(param)) {
            stop ('Only1Allele must be either TRUE or FALSE',call.=FALSE)
        }
    }

    if (name == "simul_output") {
        tmpdir <- createtempdir(param)
        return(tmpdir)
    }

    if (name == "simulation_script") {
        if (is.null(param)) {
            stop('simulation_script cannot be NULL')
        }

        if (!file.exists(param)) {
            stop('Cannot find <', param, '> please include complete path to the simulation script')
        }
    }

    if (name == "skipScriptRun") {
        if (!is.logical(param)) {stop("cannot continue, skipScriptRun has to be either TRUE or FALSE")}
    }

    if (name == "Iter") {
        if (class(param) != 'numeric') {
            stop ('Iter must be a numeric value',call.=FALSE)
        }
    }

    if (name == "conf_level") {
        if (class(param) != 'numeric') {
            stop ('conf_level must be a numeric value',call.=FALSE)
        }
    }

    if (name == "cores") {
        if (class(param) != 'numeric') {
            stop ('cores must be a numeric value (e.g. cores = 4)',call.=FALSE)
        }
    }

    if (name == "RMcorrection") {
        if (class(param) != 'logical') {
            stop ('RMcorrection must be a logical value',call.=FALSE)
        }
    }

    if (name == "RAFcorrection") {
        if (class(param) != 'logical') {
            stop ('RAFcorrection must be a logical value',call.=FALSE)
        }
    }
    #if (name == "RAF_tr") {
    #   if (class(param) == 'numeric' & length(param) != 2) {
    #       stop ('RAF_tr must be a numeric vector of length 2',call.=FALSE)
    #   }
    #   if (class(param) != 'numeric' & !is.null(param)) {
    #       stop ('RAF_tr must be a numeric vector of length 2',call.=FALSE)
    #   }
    #}

    if (name == "get.what") {
        opts<-c("samples", "param", "mergedCounts", "alleleCountsPerBam", "assayedVar", "biasTable")
         if(!is.character(param) || length(param)!=1 || !(param %in% opts))
            stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
    }

    if (name == "plot.what") {
        opts<-c("simulation_stats","barplot_per_group","boxplot_per_filter","overall_pie","pie","barplot","boxplot")
         if(!is.character(param) || length(param)!=1 || !(param %in% opts))
            stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
    }

}





