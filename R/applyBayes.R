#BaalChIP: functions to apply bayesian framework for allele-specific detection
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz

############### SNP_hit_Peaks: sum read counts ##################


getMergedRepData <- function(id_col, peak_id, SNP_hit_Peaks) {
    nrTargets <- length(peak_id)

    if (id_col == nrTargets) {
            id_REF <- seq(peak_id[id_col]+1,ncol(SNP_hit_Peaks), by=2)
            id_TOT <- seq(peak_id[id_col]+1,ncol(SNP_hit_Peaks), by=1)
    }else {
            id_REF <- seq(peak_id[id_col]+1,peak_id[id_col+1]-1, by=2)
            id_TOT <- seq(peak_id[id_col]+1,peak_id[id_col+1]-1, by=1)
    }

    if (length(id_REF) ==1) {
            return(cbind(SNP_hit_Peaks[,peak_id[id_col]], SNP_hit_Peaks[,id_REF], rowSums(SNP_hit_Peaks[,id_TOT], na.rm=TRUE)))
    } else {
            return(cbind(SNP_hit_Peaks[,peak_id[id_col]], rowSums(SNP_hit_Peaks[,id_REF], na.rm=TRUE), rowSums(SNP_hit_Peaks[,id_TOT], na.rm=TRUE)))
    }

}

Sum_read_counts <- function(SNP_hit_Peaks) {

    COLname <- colnames(SNP_hit_Peaks)
    if(COLname[1] !="ID") {
        stop("Please check the column names of SNP table: The first colname should be: ID ")
    }

    peak_id <- which(grepl("score", COLname))
    mergedReps <- lapply(seq_len(length(peak_id)), getMergedRepData, peak_id=peak_id, SNP_hit_Peaks=SNP_hit_Peaks)
    mergedReps <- data.frame("ID" = SNP_hit_Peaks[,1], do.call(cbind, mergedReps), stringsAsFactors=FALSE)

    COLname <- colnames(SNP_hit_Peaks)
    colnames(mergedReps) <- c("ID", unlist(lapply(COLname[peak_id], function (x) {c(x,"Ref","Total")})))

    return(mergedReps)
}


############## apply Bayes #################

applyBayes <- function(snp_start, snp_end, Iter, TF_num,SNP_hit_Peaks_sum, SNP_Bias) {

    ################################################
    # beta binomial functions for likelihood
    sigmoid <- function(x) {1/(1 + exp(-x))}

    invsigmoid <- function(z) {log(z) - log(1-z)}

    sigmoidln <- function(x) {-log(1.0 + exp(-x))}

    bound_sigmoid <- function(x) {(1-1e-10)*(sigmoid(x)-0.5) + 0.5}

    normpdfln <- function(x,m,std) {sum(-0.5*log(2*pi) - log(std) - 0.5*((x-m)/std)**2)}

    logit <- function(x) {log(x) - log(1-x)}

    sigmoid_trans <- function(x, opt) {
        if (opt == "tran_x") {logit(x)}
        else if (opt == "orig_x") {bound_sigmoid(x)}
        else if (opt == "det_jacob") {log(x) + log(1-x)}
        else if (opt == "tran_grad") {x*(1-x)}
        else if (opt == "jacob_grad") {1-2*x}
    }


    log_factorial <- function(n) {
        lgamma(n+1)
    }

    log_binomial_coefficient <- function(n,x) {
        log_factorial(n) - log_factorial(x) - log_factorial(n - x)
    }

    log_beta <- function(a,b) {
        if (lgamma(a) == Inf) {
        cat("a=",a)
        }else if (lgamma(b) == Inf) {
        cat("b=",b)
        }
        lgamma(a) + lgamma(b) - lgamma(a + b)
    }

    log_beta_binomial_pdf2 <- function(X, a, b){
        if (class(X) == "list") {X <- unlist(X)}
        if (X[1] == 1) {
            x = X[2]
            n = X[3]
            return(log(pi) + log_binomial_coefficient(n, x) + log_beta(a + x, b + n - x) - log_beta(a, b))
        }else {return(0)}
    }

    log_beta_binomial_pdf <- function(x, n, a, b){
        log_binomial_coefficient(n, x) + log_beta(a + x, b + n - x) - log_beta(a, b)
    }
    ################################################
    log_genotype_llh2 <- function(X_count, bias_in, allele_bias, precision){

        mu <- bias_in
        xi <- (allele_bias * mu)/(1-allele_bias-mu + 2*allele_bias*mu)


        param_a <- xi * precision
        param_b <- (1 - xi) * precision

        X_count <- lapply(seq(1, length(X_count), by=3), function(x){X_count[x: (x+2)]})
        r <- unlist(lapply(X_count, log_beta_binomial_pdf2, a=param_a, b=param_b))
        return(r)
    }

    log_genotype_llh <- function(A_count, total_count, bias_in, allele_bias, precision){

        mu <- bias_in
        xi <- (allele_bias * mu)/(1-allele_bias-mu + 2*allele_bias*mu)


        param_a <- xi * precision
        param_b <- (1 - xi) * precision
        return(log_beta_binomial_pdf(A_count, total_count, param_a, param_b))
    }

    ################################################
    MH_iter <- function(Iter,TF_num,SNP_hit_Peaks_sum, SNP_Bias, SNP_id) {
        if (identical(SNP_Bias[SNP_id,"RAF"],0)) { SNP_Bias[SNP_id,"RAF"] <- 0.01}
        if (identical(SNP_Bias[SNP_id,"RAF"],1))  { SNP_Bias[SNP_id,"RAF"] <- 0.99}

        bias <- matrix(0,Iter,1)
        bias[1] <- 0.5
        llh <- matrix(0,Iter,1)
        llh[1] <- log_pro_density_bias(bias[1],TF_num,SNP_hit_Peaks_sum, SNP_Bias, SNP_id)

        for (iter in 2:Iter) {
            set.seed(iter)
            bias_new <- bias[iter-1] + rnorm(1,mean=0, sd = 0.2)
            llh_new <- log_pro_density_bias(bias_new,TF_num,SNP_hit_Peaks_sum, SNP_Bias, SNP_id)
            ratio <- llh_new -llh[iter-1]

            U <- log(runif(1))


            if(U < min(0,ratio)) {
                bias[iter] <- bias_new
                llh[iter] <- llh_new
            }else{
                bias[iter] <- bias[iter-1]
                llh[iter] <- llh[iter-1]
            }
        }

    return(bias)
    }
    ################################################
    log_pro_density_bias <- function(allele_bias,TF_num,SNP_hit_Peaks_sum, SNP_Bias, SNP_id) {
        temp <- matrix(0,1,TF_num)
        i <- SNP_id
        precision <- 1000
        pi <- 0.5
        bias_in <- SNP_Bias[i,"RAF"]

        for (id_tf in c(1:TF_num) ) {

          if (SNP_hit_Peaks_sum[i, 2+3*(id_tf-1)]==1) {
            Ref_count <- SNP_hit_Peaks_sum[i,2+3*(id_tf-1)+1]
            total_count <- SNP_hit_Peaks_sum[i, 2+3*(id_tf-1)+2]
            temp[id_tf] <- log(pi) + log_genotype_llh(Ref_count, total_count, bias_in,
                                                     allele_bias, precision) # likelihood for each TF
          }
        }

        # refBias as model prior
        mu <- SNP_Bias[i,"RMbias"]
        var <- 0.05
        alpha <- ( (1-mu)/var- 1/mu) * mu^2
        beta <- alpha*(1/mu - 1)
        prior <- log(dbeta(allele_bias,alpha,beta))

        # total posterior = prior * likelihood
        total_pos <- sum(temp) + prior + log(bias_in/(allele_bias*bias_in+(1-allele_bias)*(1-bias_in)) - allele_bias*bias_in*(2*bias_in-1)/(allele_bias*bias_in+(1-allele_bias)*(1-bias_in))^2)
        return(total_pos)

   }

   log_pro_density_bias2 <- function(allele_bias,TF_num,SNP_hit_Peaks_sum, SNP_Bias, SNP_id) {
        i <- SNP_id
        precision <- 1000
        pi <- 0.5
        bias_in <- SNP_Bias[i,"RAF"]

        X_count <- SNP_hit_Peaks_sum[i, 2:ncol(SNP_hit_Peaks_sum)]
        temp <- log_genotype_llh2(X_count,bias_in, allele_bias, precision)
        temp <- matrix(data=unlist(temp), nrow=1, ncol=TF_num)

        # refBias as model prior
        mu <- SNP_Bias[i,"RMbias"]
        var <- 0.05
        alpha <- ( (1-mu)/var- 1/mu) * mu^2
        beta <- alpha*(1/mu - 1)
        prior <- log(dbeta(allele_bias,alpha,beta))

        # total posterior = prior * likelihood
        total_pos <- sum(temp) + prior + log(bias_in/(allele_bias*bias_in+(1-allele_bias)*(1-bias_in)) - allele_bias*bias_in*(2*bias_in-1)/(allele_bias*bias_in+(1-allele_bias)*(1-bias_in))^2)
        return(total_pos)

   }

  ############## parellel computing #################
  foreach(SNP_id=snp_start:snp_end, .combine='cbind')%dopar%
        MH_iter(Iter,TF_num,SNP_hit_Peaks_sum, SNP_Bias,SNP_id)
}


runBayes <- function(counts, bias, Iter=5000, conf_level=0.99, cores=4)
{
    ### RunBayes for each cell/individual
    #suppressPackageStartupMessages(require(doMPI))
    #cl <- startMPIcluster(count=10)
    #registerDoMPI(cl)

    #suppressPackageStartupMessages(require("foreach"))
    #suppressPackageStartupMessages(require(doParallel))
    #suppressPackageStartupMessages(require(coda))
    cores = cores# how many cores to use parallelly
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    ##------calculate args
    TF_num = sum(grepl("score", colnames(counts))) # number of TFs
    START = 1   # start SNP
    END = nrow(counts) # end SNP. Using the whole table:

    ##------ pool data between replicates
    counts.pooled <- Sum_read_counts(counts)

    ##------run bayesian model
    print(system.time(
    iter_matrix <- applyBayes(START,END,Iter,TF_num,counts.pooled,bias)
    ))

    ##------generate report
    # significant level
    threshold_lower = 0.4
    threshold_upper = 0.6

    # parameters for MCMC plot
    burnin = 0.2*Iter
    maxlag = 150
    SNP_check   = 4
    Bayes_report <- Bayesian_report(iter_matrix,conf_level,threshold_lower,threshold_upper,burnin,maxlag,SNP_check, counts)

    stopCluster(cl)
    #mpi.quit()
    return(Bayes_report)
}
