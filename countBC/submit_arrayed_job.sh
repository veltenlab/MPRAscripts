#!/bin/bash
#$ -N count
#$ -cwd 
#$ -e count.$TASK_ID.err
#$ -o count.$TASK_ID.out
#$ -q short-sl7
#$ -t 1-4
#$ -l virtual_free=50G,h_rt=5:00:00

#very memory intense because it remembers all the many UMIs for each barcode...
#if this does not work i have to srt the fastq files on hard disk by barcode seq :-/

module use /software/as/el7.2/EasyBuild/CRG/modules/all

#Nextseq2000 results and all CRS-BC associations from the NextSeq association run
FASTQDIR=/users/lvelten/sequencing_data/Chelsea_Szu-Tu/2022-08-31_AAANW5HHV/VELTENLAR_44
ASSIGNMENT="/users/lvelten/lvelten/Analysis/SCG4SYN/LibA/CRS_BC_ass/myPipeline/outNextseq/mapped.csv"
OUT="/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/out/EoBaso/"

#Sample list
INDEXX=$(echo "${SGE_TASK_ID} - 1" | bc )
SAMPLES=("3EDNA1_lib_10156AAD_TGAACTTCTT" "3EDNA2_lib_10157AAD_GAAGCTGAAG" "3ERNA1_lib_10158AAD_GATCCGGTTG" "3ERNA2_lib_10159AAD_CGTCCGGCAT")
EXP=$(echo ${SAMPLES[${INDEXX}]} | perl -pe 's/_.+//')
DNA1=$(echo ${SAMPLES[${INDEXX}]})
#Lane info
lane="001.fastq.gz"

#Specify each sample
FWD="$FASTQDIR/$DNA1""_R1_""$lane";
REV="$FASTQDIR/$DNA1""_R3_""$lane";
UMI="$FASTQDIR/$DNA1""_R2_""$lane";
echo $EXP
echo $DNA1
perl count_BC.pl $UMI $FWD $REV $ASSIGNMENT $OUT/$EXP.csv.gz $OUT/$EXP.umistat.txt $OUT/$EXP.nonmapped.csv.gz 2> $OUT/$EXP.log
