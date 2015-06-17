#BaalChIP: functions for argument checking
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

readsamplesheet <- function(samplesheet, .CHECKS=TRUE) {
	# Read samplesheet
    try(samples <- read.delim(samplesheet,stringsAsFactors=F), silent = TRUE)
    if (!exists("samples")) { stop('could not read samplesheet',call.=FALSE) } 
    
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
        for (rownr in 1:nrow(samples)) {
            x <- samples[rownr,]
            if (!file.exists(x[["bam_name"]]))
                {stop(paste('BAM file does not exist:', x[["bam_name"]]),call.=FALSE)}
            if (!file.exists(paste0(x[["bam_name"]],".bai")))
            {stop(paste('BAM index file does not exist:', paste0(x[["bam_name"]],".bai")),call.=FALSE)}
            if (!file.exists(x[["bed_name"]]))
                {stop(paste('BED file does not exist:', x[["bed_name"]]),call.=FALSE)}
        }
    }    
    
    #return
    cat("-samplesheet checks: OK!\n")
    return(samples)
}

readhettables <- function(hets, .CHECKS=TRUE) {
    if (.CHECKS) {
        for (filename in hets) {
            if (!file.exists(filename)) { stop(paste('hetSNP file does not exist:',filename),call.=FALSE) } 
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

checkmatchingnames <- function(names1, names2) {
	if (length(names1) > length(names2)) {warning(paste("failed for", setdiff(names1, names2)),call.=FALSE)}
    if (length(names2) > length(names1)) {stop(paste("failed for", setdiff(names2, names1)),call.=FALSE)}
    if (!all(names2 %in% names1)) {stop(paste("failed for", names2, names1),call.=FALSE)}
    
}

BaalChIP.checks <- function(name, param, .CHECKS= TRUE){
	
	if (name == "samplesheet")
		{
		samples <- readsamplesheet(param, .CHECKS = .CHECKS)
		return(samples)
	}
	
	if (name == "hets")
		{
		readhettables(param, .CHECKS = .CHECKS)
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
		trycreatedir(param)
	}
	
	if (name == "Iter") {
		if (class(param) != 'numeric') {
			stop ('Iter must be a numberic value',call.=FALSE)
		}
	}
	
	if (name == "conf_level") {
		if (class(param) != 'numeric') {
			stop ('conf_level must be a numberic value',call.=FALSE)
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
	#	if (class(param) == 'numeric' & length(param) != 2) {
	#		stop ('RAF_tr must be a numeric vector of length 2',call.=FALSE)
	#	}
	#	if (class(param) != 'numeric' & !is.null(param)) {
	#		stop ('RAF_tr must be a numeric vector of length 2',call.=FALSE)
	#	}
	#}
	
	if (name == "get.what") {
		opts<-c("samples", "param", "mergedCounts", "alleleCountsPerBam")
         if(!is.character(param) || length(param)!=1 || !(param %in% opts))
      		stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
	}
	
	if (name == "plot.what") {
		opts<-c("simulation_stats","barplot_per_group","boxplot_per_filter","overall_pie","pie","barplot","boxplot")
         if(!is.character(param) || length(param)!=1 || !(param %in% opts))
      		stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
	}
	
}
	




