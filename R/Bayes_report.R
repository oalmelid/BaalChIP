# BaalChIP: bayesian_report function to report the bayesian analysis result per SNP
# Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

boundary_titles <- function(threshold) {
    c(
        gettextf("conf_%.2f_lower", threshold),
        gettextf("conf_%.2f_upper", threshold)
    )
}

bayesian_report <- function(SNP_id,
                            iter_matrix,
                            conf_level,
                            central_threshold,
                            burnin,
                            maxlag) {
    traces <- iter_matrix
    mcmc_traces <- as.mcmc(traces[-(1:burnin), ])

    conf_itvals <- unlist(lapply(conf_level, function(x) {
        HPDinterval(mcmc_traces, prob = x)
    }))
    num_columns <- 3 + length(conf_itvals)

    conf_itval <- matrix(NA, 1, num_columns)
    conf_itval[1, 1:2] <- conf_itvals[1:2]
    conf_itval[1, 3] <- (0.5 + central_threshold <= conf_itval[1])
    conf_itval[1, 4] <- (0.5 - central_threshold >= conf_itval[2])
    stat_summaries <- summary(mcmc_traces)[["statistics"]]
    conf_itval[1, 5] <- stat_summaries["SD"]

    conf_columns <- c()
    if (length(conf_level) > 1) {
        conf_itval[1, 6:num_columns] <- conf_itvals[3:length(conf_itvals)]
        conf_columns <- lapply(
            conf_level[2:length(conf_level)],
            boundary_titles
        )
    }

    conf_itval_df <- data.frame(c(SNP_id), conf_itval, stringsAsFactors = FALSE)

    colnames(conf_itval_df) <- c(
        "ID",
        "Bayes_lower",
        "Bayes_upper",
        "Bayes_sig_A",
        "Bayes_sig_B",
        "Bayes_SD",
        unlist(conf_columns)
    )

    return(conf_itval_df)
}
