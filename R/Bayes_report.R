#BaalChIP: Bayesian_report function to report the bayesian analysis result
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


Bayesian_report <- function(SNP_id, iter_matrix,conf_level,threshold_lower,threshold_upper,burnin,maxlag,SNP_check){
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
    
    conf_itval <- matrix(NA,1,5)
    conf_itval[SNP_id,1:2] = HPDinterval(mcmc_traces[[SNP_id]], prob = conf_level)
    conf_itval[SNP_id,3] = (threshold_upper <= conf_itval[SNP_id,1])
    conf_itval[SNP_id,4] = (threshold_lower >= conf_itval[SNP_id,2])
    stat_summaries <- summary(mcmc_traces[[SNP_id]])[["statistics"]]
    conf_itval[SNP_id,5] = stat_summaries["SD"]
    
    conf_itval_df <- data.frame(c(SNP_id), conf_itval, stringsAsFactors=FALSE)
    colnames(conf_itval_df) <- c("ID","Bayes_lower", "Bayes_upper","Bayes_sig_A", "Bayes_sig_B", "Bayes_SD")
    
    return(conf_itval_df)
}
