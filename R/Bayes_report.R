#BaalChIP: Bayesian_report function to report the bayesian analysis result
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


Bayesian_report <- function(iter_matrix,conf_level,threshold_lower,threshold_upper,burnin,maxlag,SNP_check,SNP){
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

    conf_itval[SNP,1:2] = HPDinterval(mcmc_traces[[SNP]], prob = conf_level)
    cconf_itval[SNP,3] = (threshold_upper <= conf_itval[SNP,1])
    cconf_itval[SNP,4] = (threshold_lower >= conf_itval[SNP,2])
    stat_summaries <- summary(mcmc_traces[[SNP]])[["statistics"]]
    conf_itval[SNP,5] = stat_summaries["SD"]
    
    conf_itval <- data.frame(SNP, conf_itval, stringsAsFactors=FALSE)

    ##################
    colnames(conf_itval) <- c("ID","Bayes_lower", "Bayes_upper","Bayes_sig_A", "Bayes_sig_B", "Bayes_SD")
    return(conf_itval)
}
