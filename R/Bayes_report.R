#BaalChIP: Bayesian_report function to report the bayesian analysis result
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


Bayesian_report <- function(iter_matrix,conf_level,threshold_lower,threshold_upper,burnin,maxlag,SNP_check,SNP_hit_Peaks){
    #suppressPackageStartupMessages(require('coda'))
    traces <- iter_matrix

    mcmc_converter <- function(traces, burnin){
    numsamp = dim(traces)[1]
    numvar = dim(traces)[2]
    output = lapply(1:numvar, function(nvar)
      as.mcmc(traces[(burnin+1):numsamp, nvar]) )
    return( output )
    }


    mcmc_traces = mcmc_converter(traces, burnin)
    #autocorr.plot(mcmc_traces[[nvar]], maxlag)
    #traceplot(mcmc_traces[[nvar]])
    #densplot(mcmc_traces[[nvar]])
    #geweke.diag(mcmc_traces[[nvar]])

    conf_itval <- matrix(NA,ncol(traces),4)
    for (SNP in 1:ncol(traces)) {
    conf_itval[SNP,1:2] = HPDinterval(mcmc_traces[[SNP]], prob = conf_level)
    conf_itval[SNP,3] = (threshold_upper <= conf_itval[SNP,1])
    conf_itval[SNP,4] = (threshold_lower >= conf_itval[SNP,2])

    }
    conf_itval <- data.frame(as.character(SNP_hit_Peaks[,1]),conf_itval, stringsAsFactors=FALSE)

    ##################
    colnames(conf_itval) <- c("ID","Bayes_lower", "Bayes_upper","Bayes_sig_A", "Bayes_sig_B")
    return(conf_itval)
    #save(Bayes_binom, file=paste("results/AAF_refbias_new_post/bayes_binom/",cell_line,"_table.RData",sep=""))
    #write.table(Bayes_binom, file = paste("results/AAF_refbias_new_post/bayes_binom/",cell_line,"_table.csv", sep=""), quote=FALSE, sep=",",row.names=FALSE)

}
