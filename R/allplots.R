#BaalChIP: plot funtions
#Ines de Santiago, Wei Liu, Ke Yuan, Florian Markowetz


plot.simul <- function(simulation_stats, plot=TRUE) {
    #suppressPackageStartupMessages(require(ggplot2))
    p <- ggplot(simulation_stats, aes(x = readslen, y = perc_right)) +
                    geom_point(col="red", size=5) + geom_line() +
                    ylab("Right calls (%)") + xlab("read length")
    if (plot) {plot(p)} else{return(p)}
}


plot.filt.barplot <- function(filtering_stats, col=NULL, X_ORDER=NULL, addlegend=TRUE, plot=TRUE) {
    #plot filtering_stats
    if (addlegend==TRUE) {addlegend = "legend"}
    #suppressPackageStartupMessages(require(reshape2))
    #suppressPackageStartupMessages(require(ggplot2))
    filtering_stats$cellname <- rownames(filtering_stats)
    data2plot_stats <- melt(filtering_stats, id="cellname")

    if (!is.null(X_ORDER)) {data2plot_stats$cellname <- factor(data2plot_stats$cellname,levels=X_ORDER)}
    p <- ggplot(data=data2plot_stats, aes(x=cellname, y=value, fill=variable)) +
                geom_bar(stat="identity", position="fill", colour="black") +
                ylab("Filtered/Total (%)") + xlab("cell line") +
                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text=element_text(size=14), axis.title=element_text(size=14))

    if (!is.null(col)) {p <- p + scale_fill_manual(values=col, guide=addlegend)}
    if (plot) {plot(p); return(NULL)} else{return(p)}
}

plot.filt.boxplot <- function(filtering_stats, col=NULL, COLORGROUPS=NULL, addlegend=TRUE, plot=TRUE) {
    #suppressPackageStartupMessages(require(reshape2))
    #suppressPackageStartupMessages(require(doBy))
    #suppressPackageStartupMessages(require(ggplot2))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total), stringsAsFactors=FALSE)
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_perc <-  melt(filtering_perc, id="cellname")


    #plot filtering_stats
    if (length(unique(rownames(filtering_stats)))==1) {
                #warning("only 1 cell line, not able to plot properly...");
                filtering_stats <- rbind(filtering_stats,filtering_stats) #hack!
                }
    if (!is.null(COLORGROUPS)) {data2plot_perc$coltype <- COLORGROUPS[data2plot_perc$cellname]}
    if (is.null(COLORGROUPS)) {data2plot_perc$coltype <- 0; addlegend=FALSE}


    p <- ggplot(data2plot_perc, aes(x = variable, y = value)) + geom_boxplot(lwd=1.05, outlier.colour=NA)  +
                ylab("Filtered/Total (%)") + xlab("Filter") +
                geom_jitter(aes(colour = coltype),
                position = position_jitter(width = .1)) +
                theme(axis.text=element_text(size=14, angle=90, hjust=1), axis.title=element_text(size=14))

    if (addlegend==FALSE) {p <- p+guides(colour=FALSE)}

    if (plot) {plot(p); return(NULL)} else{return(p)}
}


plot.filt.pie <- function(filtering_stats, col=NULL, addlegend=TRUE, plot=TRUE) {
    #plot filtering_stats
    #suppressPackageStartupMessages(require(ggplot2))
    #suppressPackageStartupMessages(require(reshape2))
    #suppressPackageStartupMessages(require(doBy))

    if (addlegend==TRUE) {addlegend = "legend"}
    if (length(unique(rownames(filtering_stats)))==1) {
                #warning("only 1 cell line, piechart with total counts");
                filtering_stats <- rbind(filtering_stats,filtering_stats) #hack!
                }

    #suppressPackageStartupMessages(require(reshape2))
    #suppressPackageStartupMessages(require(doBy))
    #suppressPackageStartupMessages(require(ggplot2))
    total <- rowSums(filtering_stats)
    filtering_perc <- data.frame(apply(filtering_stats,2, function (x) 100 * x / total), stringsAsFactors=FALSE)
    filtering_perc$cellname <- rownames(filtering_perc)
    data2plot_perc <-  melt(filtering_perc, id="cellname")
    meansPERC <- summaryBy(value ~ variable , data2plot_perc, FUN=c(mean))

    p <- ggplot(meansPERC, aes(x=factor(0),fill=factor(variable),weight=value.mean)) +
                        geom_bar(width=1) +
                        coord_polar(theta='y')

    if (!is.null(col)) {p <- p + scale_fill_manual(values=col, guide=addlegend)}
    if (plot) {plot(p); return(NULL)} else{return(p)}
}

plotfilters <- function(stats, what=c("barplot_per_group","boxplot_per_filter","overall_pie"), col=NULL, X_ORDER=NULL, COLORGROUPS=NULL, addlegend=TRUE, plot=TRUE){
    cbPalette1 <- c("#000000", "#E69F00", "#69DA12", "#009E73", "#F0E442", "#0072B2")  #Color
    cbPalette2 <- c("#000000", "#E69F00", "#69DA12", "#CC79A7", "#009E73", "#F0E442", "#0072B2")  #Color
    if (is.null(col)) {
            N <- ncol(stats$filtering_stats)
            if (N < 6) {col <- cbPalette1[c(1:(N-1),6)]}
            if (N == 6) {col <- cbPalette1}
            if (N == 7) {col <- cbPalette2}
    }
    what <- match.arg(what, c("barplot_per_group","boxplot_per_filter","overall_pie"))
    switch(what,
           barplot_per_group={
              p <- plot.filt.barplot (stats$filtering_stats, col, X_ORDER=X_ORDER, addlegend=addlegend, plot=plot)
           },
           boxplot_per_filter={
             p <- plot.filt.boxplot (stats$filtering_stats, col, COLORGROUPS=COLORGROUPS, addlegend=addlegend, plot=plot)
           },
           overall_pie={
             p <- plot.filt.pie     (stats$filtering_stats, col, addlegend=addlegend, plot=plot)
           })

    if (!plot) {return(p)}else{return(NULL)}
}

plotadjustment <- function(report, col=c( "green3","gray50")) {

    #require(reshape2)
    #require(ggplot2)
    #require(scales)

    col1 <- alpha(col, .7)
    col2 <- alpha(col, .5)

    group_names <- unique(names(report))
    for(group_name in group_names) {report[[group_name]][["group_name"]] <- group_name}

    data2plot <- data.frame(do.call("rbind", report), stringsAsFactors=FALSE)
    rownames(data2plot) <- NULL
    data2plot$group_name <- factor(data2plot$group_name, levels=group_names)
    data2plot2 <- data2plot[,c("group_name","ID","AR","Corrected.AR"), drop=FALSE]
    data2plot2 <- melt(data2plot2, id=c("ID","group_name"))

    a <- ggplot(data=data2plot2, aes(x=value, fill=variable, colour=variable)) +
        geom_density(adjust=1.5, alpha=0.2) +
        scale_fill_manual(values = col1) +
        scale_colour_manual(values = col2) +
        facet_wrap(~group_name, scales="free", ncol=8) +
        theme(axis.text = element_text(color="gray25", size=6))
  plot(a)
}
