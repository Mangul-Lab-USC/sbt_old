#!/bin/bash


if [ $# -lt 2 ]
    then
    echo "********************************************************************"
    echo "Seeing Beyond the Target (SBT) - computational protocol to process off target reads in WXS experiments"
    echo "This  was written by Serghei Mangul"
    echo "********************************************************************"
    echo "1 <bam>   - bam file with mapped and unmapped reads"
    echo "2 <outdir>  - dir to save the results"
    echo "--------------------------------------"
    exit 1
    fi

. /u/local/Modules/default/init/modules.sh
module load bowtie2
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
SAMPLE=$(basename $fastq | awk -F ".bam" '{print $1}')
DIR_CODE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
bam=$1


#1. Extract unmapped reads
echo "Extract unmapped reads from " $1
$samtools view -f 0x4 -bh $bam >${SAMPLE}.unmapped.fastq

#2. Extract reads from GL000220.1 
$samtools view -bh ${bam} GL000220.1 | $samtools fastq - >${SAMPLE}.GL000220.fastq


#3 5S chr1:226743523–231,781,906???  last coordinate in one bam file249,240,475
python ${DIR_CODE}/5S.py $bam ${SAMPLE}.5S.fastq


# 4. merge unmapped + GL000220.1 + 5S and map to 45S and 5S
module load bcftools

cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.GL000220.fastq ${SAMPLE}.5S.fastq >${SAMPLE}.cat.rDNA.fastq


bowtie2  -x ${DIR_CODE}/../db.human/rDNA --end-to-end -D 15 -R 2 -L 22 -i S,1,1.15 ${SAMPLE}.cat.rDNA.fastq | $samtools view -F 4 -bh - | $samtools sort - >${sample}.sort.rDNA.bam

$samtools index ${sample}.sort.rDNA.bam
$samtools mpileup -uf  ${DIR_CODE}/../db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.bam >${SAMPLE}.rDNA.pileup



# --------------------------------------------------------
# 5.merge MT and unmapped and map to MT

$samtools view -bh $bam MT | $samtools fastq - >${SAMPLE}.MT.fastq
cat ${SAMPLE}.MT.fastq ${SAMPLE}.unmapped.fastq >${SAMPLE}.cat.MT.fastq

bowtie2  -x ${DIR}/../db.human/MT --end-to-end -D 15 -R 2 -L 22 -i S,1,1.15 ${SAMPLE}.cat.MT.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.MT.bam

$samtools index ${SAMPLE}.sort.MT.bam


$samtools mpileup -uf ${DIR_CODE}/../db.human/MT.fasta ${SAMPLE}.sort.MT.bam >${SAMPLE}.MT.pileup

#bcftools  call -mv -Oz >${SAMPLE}.MT.bcf


# --------------------------------------------------------
#6 extrcat TCR/BCR reads run imrep
#IGH
$samtools view $bam 14:106032614-107288051 |  samtools fastq - >${SAMPLE}.ireceptor.fastq
#IGK
$samtools view $bam 2:89156874-89630436 |  samtools fastq - >>${SAMPLE}.ireceptor.fastq
#IGL
$samtools view $bam 22:22380474-23265085 |  samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRA
$samtools view $bam 14:22090057-23021075 |  samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRB
$samtools view $bam 7:141998851-141998851 |  samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRG
$samtools view $bam 7:38279625-38407656 |  samtools fastq - >>${SAMPLE}.ireceptor.fastq

cat ${SAMPLE}.ireceptor.fastq ${SAMPLE}.unmapped.fastq >${SAMPLE}.cat.ireceptor.fastq

bowtie2  -x ${DIR}/../db.human/IMGT/receptor.genes --end-to-end -D 15 -R 2 -L 22 -i S,1,1.15 ${SAMPLE}.cat.ireceptor.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.ireceptor.bam

$samtools index ${SAMPLE}.sort.ireceptor.bam

$samtools mpileup -uf ${DIR}/../db.human/IMGT/receptor.genes.fasta ${SAMPLE}.sort.ireceptor.bam>${SAMPLE}.ireceptor.pileup



#imrep







  
