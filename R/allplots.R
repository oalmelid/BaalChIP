#BaalChIP plot funtions

#plot.simul <- function(simulation_stats) {
#    suppressPackageStartupMessages(require(ggplot2))
#    p <- ggplot(simulation_stats, aes(x = readslen, y = perc_right)) + 
#                    geom_point(col="red", size=5) + geom_line() + ylab("Right calls (%)") + xlab("read length") + ylim(c(40,100)) 
#    plot(p)    
#}

plot.filt.barplot <- function(filtering_stats, col=NULL, X_ORDER=NULL, addlegend=TRUE) {
    #plot filtering_stats
    if (addlegend==TRUE) {addlegend = "legend"}
    suppressPackageStartupMessages(require(reshape2))
    suppressPackageStartupMessages(require(ggplot2))
    filtering_stats$cellname <- rownames(filtering_stats)
    data2plot_stats <- melt(filtering_stats, id="cellname")
    
    if (!is.null(X_ORDER)) {data2plot_stats$cellname <- factor(data2plot_stats$cellname,levels=X_ORDER)}
    p <- ggplot(data=data2plot_stats, aes(x=cellname, y=value, fill=variable)) + 
                geom_bar(stat="identity", position="fill", colour="black") +
                ylab("Filtered/Total (%)") + xlab("cell line") + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                axis.text=element_text(size=14), axis.title=element_text(size=14))
    
    if (!is.null(col)) {p <- p + scale_fill_manual(values=col, guide=addlegend)}
    plot(p)
}

plot.filt.boxplot <- function(filtering_stats, col=NULL, COLORGROUPS=NULL, addlegend=TRUE) {
    suppressPackageStartupMessages(require(reshape2))
    suppressPackageStartupMessages(require(doBy))
    suppressPackageStartupMessages(require(ggplot2))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total))
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_perc <-  melt(filtering_perc, id="cellname")
    
    
    #plot filtering_stats
    if (length(unique(rownames(filtering_stats)))==1) {
                print("only 1 cell line, not able to plot properly..."); 
                filtering_stats <- rbind(filtering_stats,filtering_stats) #hack!
                }
    if (!is.null(COLORGROUPS)) {data2plot_perc$coltype <- COLORGROUPS[data2plot_perc$cellname]}
    if (is.null(COLORGROUPS)) {data2plot_perc$coltype <- 0; addlegend=FALSE}
    
    
    p <- ggplot(data2plot_perc, aes(x = variable, y = value)) + geom_boxplot(lwd=1.05,outlier.shape=NA)  +
                ylab("Filtered/Total (%)") + xlab("Filter") +
                geom_jitter(aes(colour = coltype), 
                position = position_jitter(width = .1)) +
                theme(axis.text=element_text(size=14, angle=90, hjust=1), axis.title=element_text(size=14))
                
    if (addlegend==FALSE) {p <- p+guides(colour=FALSE)}
                
    plot(p)
}

plot.filt.pie <- function(filtering_stats, col=NULL, addlegend=TRUE) {
    #plot filtering_stats
    suppressPackageStartupMessages(require(ggplot2))  
    suppressPackageStartupMessages(require(reshape2))
    suppressPackageStartupMessages(require(doBy))
    
    if (addlegend==TRUE) {addlegend = "legend"}
    if (length(unique(rownames(filtering_stats)))==1) {
                print("only 1 cell line, piechart with total counts"); 
                filtering_stats <- rbind(filtering_stats,filtering_stats) #hack!
                }
    
    suppressPackageStartupMessages(require(reshape2))
    suppressPackageStartupMessages(require(doBy))
    suppressPackageStartupMessages(require(ggplot2))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total))
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_perc <-  melt(filtering_perc, id="cellname")
    meansPERC <- summaryBy(value ~ variable , data2plot_perc, FUN=c(mean))
    
    p <- ggplot(meansPERC, aes(x=factor(0),fill=factor(variable),weight=value.mean)) +
                        geom_bar(width=1) + 
                        coord_polar(theta='y')
    
    if (!is.null(col)) {p <- p + scale_fill_manual(values=col, guide=addlegend)}
    plot(p)
}

plotfilters <- function(stats, what=c("barplot_per_group","boxplot_per_filter","overall_pie"), col=NULL, X_ORDER=NULL, COLORGROUPS=NULL, addlegend=TRUE){
    cbPalette1 <- c("#000000", "#E69F00", "#69DA12", "#009E73", "#F0E442", "#0072B2")  #Color
    cbPalette2 <- c("#000000", "#E69F00", "#69DA12", "#CC79A7", "#009E73", "#F0E442", "#0072B2")  #Color
    if (is.null(col)) {
            N <- ncol(stats$filtering_stats)
            if (N < 6) {col <- cbPalette1[c(1:(N-1),6)]}
    		if (N == 6) {col <- cbPalette1}
    		if (N == 7) {col <- cbPalette2}
    }
    if (what == "barplot_per_group")     { plot.filt.barplot (stats$filtering_stats, col, X_ORDER=X_ORDER, addlegend=addlegend) }
    if (what == "boxplot_per_filter")   { plot.filt.boxplot (stats$filtering_stats, col, COLORGROUPS=COLORGROUPS, addlegend=addlegend) }
    if (what == "overall_pie")          { plot.filt.pie     (stats$filtering_stats, col, addlegend=addlegend) }
}

