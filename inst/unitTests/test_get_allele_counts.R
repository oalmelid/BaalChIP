test_get_allele_counts <- function() {
    bamfile <- file.path(system.file("test",package="BaalChIP"), "bamFiles","MCF7_cFOS_Rep1.bam")
    snpfile <- file.path(system.file("test",package="BaalChIP"), "MCF7_hetSNP.txt")
    bedname <- file.path(system.file("test",package="BaalChIP"), "bedFiles","MCF7_cFOS.bed")
    snp.ranges <- BaalChIP:::filter_sigi(snpfile=snpfile, bedfile = bedname)

    res <- BaalChIP:::get_allele_counts(bamfile, snp.ranges)
    checkEquals(ncol(res), 10)
    checkEquals(nrow(res), 94)
    checkEquals(class(res), "data.frame")
    checkEquals(sum(res$REF.counts), 267)
    checkEquals(sum(res$ALT.counts), 317)
}
