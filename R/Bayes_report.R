#BaalChIP: Bayesian_report function to report the bayesian analysis result per SNP
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


Bayesian_report <- function(SNP_id,iter_matrix,conf_level,threshold_lower,threshold_upper,burnin,maxlag,SNP_check){
    #suppressPackageStartupMessages(require('coda'))
    traces <- iter_matrix
    mcmc_traces <- as.mcmc(traces[-(1:burnin),])
    
    conf_itval <- matrix(NA,1,5)
    conf_itval[1,1:2] = HPDinterval(mcmc_traces, prob = conf_level)
    conf_itval[1,3] = (threshold_upper <= conf_itval[1])
    conf_itval[1,4] = (threshold_lower >= conf_itval[2])
    stat_summaries <- summary(mcmc_traces)[["statistics"]]
    conf_itval[1,5] = stat_summaries["SD"]
    
    conf_itval_df <- data.frame(c(SNP_id), conf_itval, stringsAsFactors=FALSE)
    colnames(conf_itval_df) <- c("ID","Bayes_lower", "Bayes_upper","Bayes_sig_A", "Bayes_sig_B", "Bayes_SD")
    
    return(conf_itval_df)
}
