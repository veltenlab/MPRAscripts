#!/bin/bash
#$ -N map
#$ -cwd 
#$ -e mapNextseq.err
#$ -o mapNextseq.out
#$ -l virtual_free=50G,h_rt=24:00:00
#$ -q long-sl7

#when I checked the memory footprint of the perl script was <15GByte for a NextSeq MidOutput

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load BWA/0.7.17
module load SAMtools/1.10-GCC-9.3.0

FASTQDIR=/users/lvelten/sequencing_data/Chelsea_Szu-Tu/2022-02-06_AHH3FHAFX3/VELTENLAR_34

FWD="$FASTQDIR/L1GreB_15809AAC_RS_ACGTAGACCA_R1_001.fastq.gz";
REV="$FASTQDIR/L1GreB_15809AAC_RS_ACGTAGACCA_R3_001.fastq.gz";
BC="$FASTQDIR/L1GreB_15809AAC_RS_ACGTAGACCA_R2_001.fastq.gz";

OUT="/users/lvelten/lvelten/Analysis/SCG4SYN/LibA/CRS_BC_ass/myPipeline/outNextseq"

REF=/users/lvelten/lvelten/Analysis/SCG4SYN/LibA/CRS_BC_ass/data/design.fa

THREADS=1

bwa index $REF

bwa mem -t $THREADS $REF $FWD $REV | samtools view -b > $OUT/aligned.bam 2> $OUT/alignment.log

perl map_crs_bc_BAM_2.pl $BC $OUT/aligned.bam $OUT/mapped.csv 2> $OUT/mapping.log

#also need to run report.Rmd to filter the output


#this is a useful command for finding all the alignments of a given barcode
#zcat $BC | grep -B1 GTAACCACAGAGTTG | grep "M03766" | perl -pe 's/^@(\S+)\s.+/$1/' > test.txt
#samtools view $OUT/aligned.bam | grep -f test.txt 