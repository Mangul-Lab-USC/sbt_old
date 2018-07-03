#!/bin/bash



#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('inbam')
parser.add_argument('outdir')
parser.add_argument('-hg38', '--hg38', action='store_true',
                    default=False, help='Are the reads mapped to hg39 genome release. To check it please run samtools view -H <bam file> [default %(default)s]')
EOF

echo required infile: "$INFILE"
echo required outfile: "$OUTFILE"
echo $HG38


. /u/local/Modules/default/init/modules.sh
module load bowtie2
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
bam=$1
SAMPLE=$(basename $bam | awk -F ".bam" '{print $1}')

DIR=$2
mkdir $DIR

SAMPLE=$SAMPLE"/"$DIR
DIR_CODE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


#1. Extract unmapped reads
echo "---------------------------------------------------------------------------"
echo "Extract unmapped reads from " $1
$samtools view -f 0x4 -bh $bam | $samtools fastq - >${SAMPLE}.unmapped.fastq

#2. Extract reads from GL000220.1
echo "---------------------------------------------------------------------------"
echo "Extract reads from GL000220.1 and save to "${SAMPLE}.GL000220.fastq
$samtools view -bh ${bam} GL000220.1 | $samtools fastq - >${SAMPLE}.GL000220.fastq



# 4. merge unmapped + GL000220.1 + 5S and map to 45S and 5S
echo "---------------------------------------------------------------------------"
echo "Merge unmapped + GL000220.1 + 5S and map to 45S and 5S " ${SAMPLE}.cat.rDNA.fastq
module load bcftools



while read line
do
chr=$(echo $line | awk -F "," '{print $1}')
x=$(echo $line | awk -F "," '{print $2}')
y=$(echo $line | awk -F "," '{print $3}')
$samtools view -bh $bam $chr:$x-$y | $samtools fastq - >>${SAMPLE}.cat.rDNA.fastq
done<${DIR_CODE}/db.human/determine.homology/5S/homology.5S.bed


wc -l ${SAMPLE}.cat.rDNA.fastq

while read line
do
chr=$(echo $line | awk -F "," '{print $1}')
x=$(echo $line | awk -F "," '{print $2}')
y=$(echo $line | awk -F "," '{print $3}')
$samtools view -bh $bam $chr:$x-$y | $samtools fastq - >>${SAMPLE}.cat.rDNA.fastq
done<${DIR_CODE}/db.human/determine.homology/HSU13369/homology.HSU13369.bed

wc -l ${SAMPLE}.cat.rDNA.fastq

cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.GL000220.fastq ${SAMPLE}.5S.fastq >${SAMPLE}.cat.rDNA.fastq


bowtie2  -x ${DIR_CODE}/db.human/rDNA --end-to-end -D 15 -R 2 -L 22 -i S,1,1.15 ${SAMPLE}.cat.rDNA.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.rDNA.bam

$samtools index ${SAMPLE}.sort.rDNA.bam


$samtools mpileup -f  ${DIR_CODE}/db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.bam >${SAMPLE}.rDNA.pileup
$samtools mpileup -uf ${DIR_CODE}/db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.bam | bcftools  call -mv -Oz >${SAMPLE}.rDNA.bcf




# --------------------------------------------------------
# 5.merge MT and unmapped and map to MT
echo "---------------------------------------------------------------------------"
echo "merge MT and unmapped and map to MT"

$samtools view -bh $bam MT >${SAMPLE}.MT.bam
$samtools index ${SAMPLE}.MT.bam


$samtools mpileup -f ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.MT.bam >${SAMPLE}.MT.pileup
$samtools mpileup -uf ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.sort.MT.bam | bcftools  call -mv -Oz >${SAMPLE}.MT.bcf







# --------------------------------------------------------
#6 extrcat TCR/BCR reads run imrep
echo "---------------------------------------------------------------------------"
echo "extract immune reads and map to VDJ ref"
#IGH
$samtools view -bh $bam 14:106032614-107288051 |  $samtools fastq - >${SAMPLE}.ireceptor.fastq
#IGK
$samtools view -bh $bam 2:89156874-89630436 |  $samtools fastq - >>${SAMPLE}.ireceptor.fastq
#IGL
$samtools view -bh $bam 22:22380474-23265085 |  $samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRA
$samtools view -bh $bam 14:22090057-23021075 |  $samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRB
$samtools view -bh $bam 7:141998851-141998851 |  $samtools fastq - >>${SAMPLE}.ireceptor.fastq

#TCRG
$samtools view -bh $bam 7:38279625-38407656 |  $samtools fastq - >>${SAMPLE}.ireceptor.fastq

cat ${SAMPLE}.ireceptor.fastq ${SAMPLE}.unmapped.fastq >${SAMPLE}.cat.ireceptor.fastq

bowtie2  -x ${DIR_CODE}/db.human/IMGT/receptor.genes --end-to-end -D 15 -R 2 -L 22 -i S,1,1.15 ${SAMPLE}.cat.ireceptor.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.ireceptor.bam

$samtools index ${SAMPLE}.sort.ireceptor.bam

$samtools mpileup -f ${DIR_CODE}/db.human/IMGT/receptor.genes.fasta ${SAMPLE}.sort.ireceptor.bam>${SAMPLE}.ireceptor.pileup
$samtools mpileup -uf ${DIR_CODE}/db.human/IMGT/receptor.genes.fasta ${SAMPLE}.sort.ireceptor.bam | bcftools  call -mv -Oz >${SAMPLE}.ireceptor.bcf


#7 imrep
echo "---------------------------------------------------------------------------"
# 8 map to microbiome - sort of micop here! do better fungi reference TODO!
#metaphlan

# assemle contigs and map to viruses
echo "---------------------------------------------------------------------------"

echo "assemle contigs and map to viruses"
/u/home/n/nlapier2/lib/megahit_v1.1.3_LINUX_CPUONLY_x86_64-bin/megahit --k-step 10 -r ${SAMPLE}.unmapped.fastq megahit







module load bwa


bwa mem -a /u/home/s/serghei/collab/code/rop/db_human/viral.vipr/NONFLU_All.fastq ${SAMPLE}.unmapped.fastq | $samtools view -S -b -F 4 - | $samtools sort - >${SAMPLE}.virus.bam
bwa mem -a /u/home/s/serghei/collab/code/rop/db_human//fungi/fungi.ncbi.february.3.2018.fasta ${SAMPLE}.unmapped.fastq | $samtools view -S -b -F 4 - |  $samtools sort - >${SAMPLE}.fungi.bam
bwa mem -a /u/home/s/serghei/collab/code/rop/db_human/protozoa/protozoa.ncbi.february.3.2018.fasta ${SAMPLE}.unmapped.fastq | $samtools view -S -b -F 4 - | $samtools sort - >${SAMPLE}.protozoa.bam





  
python /u/home/s/serghei/collab/code/imrep/imrep.py --bam --extendedOutput ${bam} ${SAMPLE}.imrep.cdr3

