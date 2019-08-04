#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('bam')
parser.add_argument('outdir')
parser.add_argument('-hg38', '--hg38', action='store_true',default=False, help='Choose this option, if reads are mapped to hg39 genome release. To check it please run samtools view -H <bam file> [default %(default)s]')
parser.add_argument('-RNASeq', '--RNASeq', action='store_true',default=False, help=' Choose this option, if it is a RNA-Seq data[default %(default)s]')
parser.add_argument('-op', '--OncoPanel', action='store_true',default=False, help=' Choose this option, if it is a OncoPanel data[default %(default)s]')

parser.add_argument('-f', '--force', action='store_true', default=False,help='Forse [default %(default)s]')
parser.add_argument('-dev', '--dev', action='store_true', default=False,help='Keep all intermediate files [default %(default)s]')

EOF

DIR_CODE=`dirname $(readlink -f "$0")`

echo infile: "$INBAM"
echo outfile: "$OUTDIR"



#Add MiniConda to PATH if it's available.
if [ -d "$DIR_CODE/tools/MiniConda/bin" ]; then
    echo "Add MiniConda to PATH if it's available"
    export PATH="$DIR_CODE/tools/MiniConda/bin:$PATH"
fi



#Convert to absolute paths.
BAM=`readlink -m "$BAM"`
OUTDIR=`readlink -m "$OUTDIR"`
ORGANISM='human'



#Check if BAM exists.
if [ ! -e "$BAM" ]
then
    echo "Error: $BAM doesn't exist." >&2
    exit 1
fi

#Check if OUTDIR exists, then make it.
#echo $FORCE
#if [ -d "$OUTDIR" ]
#then
#    if [[ $FORCE ]]
#    then
#        rm -fr "$OUTDIR"
#    else
#        echo "Error: The directory $OUTDIR exists. Please choose a" \
#            'different directory in which to save results of the analysis, or' \
#            'use the -f option to overwrite the directory.' >&2
#        exit 1
#    fi
#fi
mkdir -p "$OUTDIR"




start=`date +%s`
echo  "Start SBT analysis ... "$start


prefix=$(basename $BAM | awk -F ".bam" '{print $1}')
PREFIX=$(basename $BAM | awk -F ".bam" '{print $1}')

SAMPLE=${OUTDIR}"/"${prefix}


megahit=${DIR_CODE}/tools/megahit/megahit
metaphlan2=${DIR_CODE}/tools/metaphlan2/metaphlan2.py
needle=${DIR_CODE}/tools/needle/needle.sh
imrep=${DIR_CODE}/tools/imrep/imrep.py


cd $OUTDIR


echo "------------- ... Extract unmapped reads from " $BAM "-------------"
samtools view -f 0x4 -bh $BAM | samtools bam2fq - >${SAMPLE}.unmapped.fastq
samtools view -bh $BAM NC_007605 | samtools fastq - > ${SAMPLE}.NC_007605.fastq
rm -fr ${SAMPLE}.NC_007605.fastq
cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.NC_007605.fastq>${SAMPLE}.cat.unmapped.fastq
rm -fr ${SAMPLE}.unmapped.fastq
UNMAPPED=${SAMPLE}.cat.unmapped.fastq


Echo "Number of unmapped reads"
wc -l  ${SAMPLE}.cat.unmapped.fastq


echo "------------- (1) Use needle to detect viruses, fungia,and protozoa -------------"
$needle --fastq $UNMAPPED ${OUTDIR}/needle.out


exit 1


echo "------------- (4) rDNA dosage -------------"
samtools view -bh ${BAM} GL000220.1 | samtools bam2fq - >${SAMPLE}.GL000220.fastq
rm -fr ${SAMPLE}.rDNA.mapped.fastq
while read line;do chr=$(echo $line | awk -F "," '{print $1}');x=$(echo $line | awk -F "," '{print $2}');y=$(echo $line | awk -F "," '{print $3}');samtools view -bh $BAM $chr:$x-$y | samtools fastq - >>${SAMPLE}.rDNA.mapped.fastq;done<${DIR_CODE}/db.human/rDNA.kmers.75.clean.filtered.bed


cat $UNMAPPED ${SAMPLE}.GL000220.fastq ${SAMPLE}.rDNA.mapped.fastq >${SAMPLE}.cat.rDNA.fastq
bowtie2  -x ${DIR_CODE}/db.human/rDNA --end-to-end ${SAMPLE}.cat.rDNA.fastq | samtools view -F 4 -bh - | samtools sort - >${SAMPLE}.sort.rDNA.bam
samtools index ${SAMPLE}.sort.rDNA.bam

samtools view -H ${SAMPLE}.sort.rDNA.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.sort.rDNA.bam | grep -v "XS:" | cat ${OUTDIR}/header.sam - | samtools view -b - > ${SAMPLE}.sort.rDNA.unique.bam
samtools depth ${SAMPLE}.sort.rDNA.unique.bam >${SAMPLE}.sort.rDNA.cov
samtools mpileup -uf ${DIR_CODE}/db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.unique.bam | bcftools  call -mv -Oz >${SAMPLE}.rDNA.bcf


rm -fr ${SAMPLE}.sort.rDNA.bam
rm -fr ${SAMPLE}.cat.rDNA.fastq ${SAMPLE}.GL000220.fastq ${SAMPLE}.rDNA.mapped.fastq ${SAMPLE}.sort.rDNA.fastq


echo "------------- (5) MT dosage and diversity -------------"
samtools view -bh $BAM MT | samtools bam2fq - >${SAMPLE}.MT.fastq
cat ${SAMPLE}.MT.fastq $UNMAPPED>${SAMPLE}.cat.MT.fastq
bowtie2  -x ${DIR_CODE}/db.human/MT --end-to-end ${SAMPLE}.cat.MT.fastq | samtools view -F 4 -bh - | samtools sort - >${SAMPLE}.sort.MT.bam

samtools index ${SAMPLE}.sort.MT.bam
samtools view -H ${SAMPLE}.sort.MT.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.sort.MT.bam | grep -v "XS:" | cat ${OUTDIR}/header.sam - | samtools view -b - > ${SAMPLE}.sort.MT.unique.bam
samtools depth ${SAMPLE}.sort.MT.unique.bam >${SAMPLE}.sort.MT.cov
samtools mpileup -uf ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.sort.MT.unique.bam | bcftools  call -mv -Oz >${SAMPLE}.MT.bcf


rm -fr ${SAMPLE}.sort.MT.bam
rm -fr ${SAMPLE}.cat.MT.fastq ${SAMPLE}.sort.MT.fastq
rm -fr ${SAMPLE}.MT.fastq



if [ $ONCOPANEL ]
then
echo "------------- (optional) TE elements. Only for OncoPanel data -------------"
samtools bam2fq $BAM | bowtie2  -x ${DIR_CODE}/db.human/repeats/repbase.fa --end-to-end - | samtools view -F 4 -bh - | samtools sort - >${SAMPLE}.sort.repeat.bam
fi


echo "------------- (6) T and B cell repetoire profiling -------------"

if [[ $HG38 ]]; then

samtools view -bh $BAM 14:105586437-106879844 |  samtools bam2fq - >${SAMPLE}.ireceptor.fastq
samtools view -bh $BAM 2:88857361-90235368 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
samtools view -bh $BAM 22:22026076-22922913 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
samtools view -bh $BAM 14:21621904-22552132 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
samtools view -bh $BAM 7:142299011-142813287 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
samtools view -bh $BAM 7:38240024-38368055 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq

else
#IGH
samtools view -bh $BAM 14:106032614-107288051 |  samtools bam2fq - >${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
#IGK
samtools view -bh $BAM 2:89156874-89630436 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
#IGL
samtools view -bh $BAM 22:22380474-23265085 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
#TCRA
samtools view -bh $BAM 14:22090057-23021075 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
#TCRB
samtools view -bh $BAM 7:141998851-141998851 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
#TCRG
samtools view -bh $BAM 7:38279625-38407656 |  samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
wc -l ${SAMPLE}.ireceptor.fastq
fi

rm -fr ${SAMPLE}.ireceptor.fastq


cat ${SAMPLE}.ireceptor.fastq $UNMAPPED >${SAMPLE}.cat.ireceptor.fastq


bowtie2  -x ${DIR_CODE}/db.human/IMGT/receptor.genes.prepared --end-to-end ${SAMPLE}.cat.ireceptor.fastq | samtools view -F 4 -bh - | samtools sort - >${SAMPLE}.sort.ireceptor.bam
samtools index ${SAMPLE}.sort.ireceptor.bam
samtools view -H ${SAMPLE}.sort.ireceptor.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.sort.ireceptor.bam | grep -v "XS:" | cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.sort.ireceptor.cov
samtools mpileup -uf ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.sort.MT.unique.bam | bcftools  call -mv -Oz >${SAMPLE}.ireceptor.bcf


#imrep
python $imrep --bam --extendedOutput ${BAM} ${SAMPLE}.imrep.cdr3




echo "------------- (7) Calculate off target coverage -------------"

while read line
do
chr=$(echo $line | awk '{print $1}')
x=$(echo $line | awk '{print $2}')
y=$(echo $line | awk '{print $3}')

n=0
n=$(samtools view -bh $BAM $chr:$x-$y | samtools  depth - | awk '{s+=$3} END {print s}')
echo $chr,$x,$y,$n >>${SAMPLE}.offtarget.cov
done<${DIR_CODE}/db.human/intergenic.regions/intergenic.regions.hg19.autosomes.bed


if [ $DEV ]
then
echo "Keep all intermediate files"
else
rm -fr ${2}/*cleaned_input.fasta
rm -fr ${2}/*fastq
rm -fr ${2}/partial*
rm -fr ${2}/*bam
rm -fr ${2}/*bai
fi

end=`date +%s`

runtime=$((end-start))

echo "It took "$runtime" to run SBT protocol on " $BAM





