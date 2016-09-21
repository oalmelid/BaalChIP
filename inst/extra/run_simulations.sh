#!/bin/sh
#BSUB -o BEdup.out.log
#BSUB -e BEdup.err.log
#BSUB -a R

#########################################################################################
# Will perform simulated set of reads for a list of SNP locations
# It is build to work within BaalChIP R package pipeline
# ---*---
#Arguments:
# ---*---
#picard: location to picard
#bowtie: location to bowtie
#genome: the location to the genome file
#genomeBychr: the location to the genome file
#simscript: complete path to get.overlaps.v2_chrXY.perl script
#readlen: read length
#snplist: a table with rsid, chr, pos, ref, alt. 
#fastaout: location of the output fasta file
#bamout: location of the output bamfile
#########################################################################################

#########################################################################################
#Enter the following arguments
#from R using 'alignmentSimulArgs' argument of filterIntbias function
#e.g.
#data(BaalObject)
#res <- filterIntbias(BaalObject,
#'       simul_output="~/simuloutput",
#'       alignmentSimulArgs=c("picard-tools-1.119","bowtie-1.1.1","genomes_test/male.hg19","genomes_test/maleByChrom"))
picard=$6
bowtie=$7
genome=$8
genomeBychr=$9
#########################################################################################

#########################################################################################
#Or 
#use alignmentSimulArgs=NULL and uncomment the following lines (using the appropriate paths):

#picard=/Volumes/groups/Research/fmlab/public_folders/InesdeSantiago/picard-tools-1.119
#bowtie=/Volumes/groups/Research/fmlab/public_folders/InesdeSantiago/bowtie-1.1.1
#genome=/Users/santia01/Desktop/genomes_test/male.hg19
#genomeBychr=/Users/santia01/Desktop/genomes_test/maleByChrom
#########################################################################################


#########################################################################################
#-----------------#
#DO NOT CHANGE
#these next five arguments are arguments of get.overlaps.v2_chrXY.perl 
#they are computed directly in R and passed to this script as arguments
#hardcoded within the R command 
simscript=$1  #simscript=/Users/santia01/Dropbox/FromHome/baal_package/BaalChIP/inst/extra/get.overlaps.v2_chrXY.perl
readlen=$2
snplist=$3
fastaout=$4
bamout=$5
#########################################################################################


#Run
readstring=`printf '%*s' "$readlen" | tr ' ' "B"`
perl $simscript $snplist $genomeBychr $readlen $readstring > $fastaout
$bowtie/bowtie --chunkmbs 512 -t --best -q -S -p 8 -l 32 -e 80 -n 2 -m 1 $genome $fastaout $fastaout.sam
java -jar -Xmx3g $picard/SortSam.jar INPUT=$fastaout.sam OUTPUT=$bamout SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
rm $fastaout.sam
mv ${bamout%.*}.bai ${bamout%.*}.bam.bai


#--chunkmbs     The number of megabytes of memory a given thread is given to store path descriptors in --best mode.
#-t             Print the amount of wall-clock time taken by each phase.
#--best         Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum (i.e. number of mismatches, or mismatches in the seed in the case of -n mode) 
#-q             With this option bowtie-build will print only error messages
#-S             Sam mode
#-p             Launch <int> parallel search threads
#-l             The "seed length"; i.e., the number of bases on the high-quality end of the read to which the -n ceiling applies.
#-e             Maximum permitted total of quality values at all mismatched read positions throughout the entire alignment, not just in the "seed". 
#-n             Maximum number of mismatches permitted in the "seed", i.e. the first L base pairs of the read 
#-m             Suppress all alignments for a particular read or pair if more than <int> reportable alignments exist for it. Reportable alignments are those that would be reported given the -n, -v, -l, -e, -k, -a, --best, and --strata options.
