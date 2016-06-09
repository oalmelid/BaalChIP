#!/bin/sh
#BSUB -o BEdup.out.log
#BSUB -e BEdup.err.log
#BSUB -a R

# Will perform simulated set of reads for a list of SNP locations
# It is build to work within BaalChIP R package pipeline
# Depends on several external scripts:
#----- 1. get.overlaps.v2_chrXY.perl (originally published by Degner et al., 2009: Bioinformatics. 2009 Dec 15;25(24):3207-12. doi: 10.1093/bioinformatics/btp579)
#----- 2. picard
#----- 3. bowtie

#Arguments to this script
#snp_dataframe: a table with rsid, chr, pos, ref, alt. 
#readlen: read length
#fastaout: location of the output fasta file
#bamout: location of the output bamfile
#location_getoverlaps: location to get.overlaps.v2_chrXY.perl script
#location_picard: location to picard
#location_bowtie: location to bowtie
#location_genome: the location to the genome file

#simscript=/lustre/fmlab/santia01/tools/Degner2009/simulate.reads/get.overlaps.v2_chrXY.perl
#picard=/lustre/fmlab/santia01/tools/picard-tools-1.47
#bowtie=/lustre/fmlab/santia01/tools/bowtie-1.1.0/bowtie
#genome=/lustre/fmlab/santia01/tools/genomes/hg19_encodeDCC/male.hg19
#genomeBychr=/lustre/fmlab/santia01/tools/genomes/hg19_encodeDCC/maleByChrom
simscript=/Users/santia01/Dropbox/FromHome/baal_package/BaalChIP/inst/extra/get.overlaps.v2_chrXY.perl
picardSortSam=/Volumes/groups/Research/fmlab/public_folders/InesdeSantiago/picard-tools-1.119/SortSam.jar
bowtie=/Volumes/groups/Research/fmlab/public_folders/InesdeSantiago/bowtie-1.1.1/bowtie
genome=/Users/santia01/Desktop/genomes_test/male.hg19
genomeBychr=/Users/santia01/Desktop/genomes_test/maleByChrom
readlen=$1
snplist=$2
fastaout=$3
bamout=$4


#Run
readstring=`printf '%*s' "$readlen" | tr ' ' "B"`
perl $simscript $snplist $genomeBychr $readlen $readstring > $fastaout
$bowtie --chunkmbs 512 -t --best -q -S -p 8 -l 32 -e 80 -n 2 -m 1 $genome $fastaout $fastaout.sam
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
