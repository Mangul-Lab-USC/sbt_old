#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('bam')
parser.add_argument('outdir')
parser.add_argument('-hg38', '--hg38', action='store_true',
                    default=False, help='Choose this option, if reads are mapped to hg39 genome release. To check it please run samtools view -H <bam file> [default %(default)s]')
parser.add_argument('-RNASeq', '--RNASeq', action='store_true',
                    default=False, help=' Choose this option, if it is a RNA-Seq data[default %(default)s]')
parser.add_argument('-f', '--force', action='store_true', default=False,
                    help='Forse [default %(default)s]')

EOF

DIR_CODE=`dirname $(readlink -f "$0")`

echo required infile: "$INBAM"
echo required outfile: "$OUTDIR"



#Add MiniConda to PATH if it's available.
if [ -d "$DIR_CODE/tools/MiniConda/bin" ]; then
    echo "Add MiniConda to PATH if it's available"
    export PATH="$DIR_CODE/tools/MiniConda/bin:$PATH"
fi



#Convert to absolute paths.
BAM=`readlink -m "$BAM"`
OUTDIR=`readlink -m "$OUTDIR"`



#Check if BAM exists.
if [ ! -e "$BAM" ]
then
    echo "Error: $BAM doesn't exist." >&2
    exit 1
fi

Check if OUTDIR exists, then make it.
echo $FORCE
if [ -d "$OUTDIR" ]
then
    if [[ $FORCE ]]
    then
        rm -fr "$OUTDIR"
    else
        echo "Error: The directory $OUTDIR exists. Please choose a" \
            'different directory in which to save results of the analysis, or' \
            'use the -f option to overwrite the directory.' >&2
        exit 1
    fi
fi
mkdir -p "$OUTDIR"




start=`date +%s`
echo  "Start SBT analysis ... "$start


prefix=$(basename $BAM | awk -F ".bam" '{print $1}')
SAMPLE=${OUTDIR}"/"${prefix}


megahit=${DIR_CODE}/tools/megahit/megahit
metaphlan2=/u/home/s/serghei/collab/code/rop/tools/metaphlan2/metaphlan2.py



echo "--------1. Extract unmapped reads from " $BAM
samtools view -f 0x4 -bh $BAM | samtools bam2fq - >${SAMPLE}.unmapped.fastq
samtools view -bh $BAM NC_007605 | samtools fastq - > ${SAMPLE}.NC_007605.fastq
rm -fr ${SAMPLE}.NC_007605.fastq
cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.NC_007605.fastq>${SAMPLE}.cat.unmapped.fastq
rm -fr ${SAMPLE}.unmapped.fastq
UNMAPPED=${SAMPLE}.cat.unmapped.fastq


echo "--------2. Use needle to detect viruses, fungia,and protozoa"
$metaphlan2 $UNMAPPED--nproc 8 --input_type fastq --bowtie2out ${SAMPLE}.metaphlan2.bowtie2.txt >${SAMPLE}.metaphlan2.txt


bwa mem -a ${DIR_CODE}/db.human/viral.vipr/NONFLU_All.fastq $UNMAPPED | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.virus.bam
bwa mem -a ${DIR_CODE}/db.human/fungi/fungi.ncbi.february.3.2018.fasta $UNMAPPED | samtools view -S -b -F 4 - |  samtools sort - >${SAMPLE}.fungi.bam
bwa mem -a ${DIR_CODE}/db.human/protozoa/protozoa.ncbi.february.3.2018.fasta $UNMAPPED | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.protozoa.bam

samtools index ${SAMPLE}.virus.bam
samtools index ${SAMPLE}.fungi.bam
samtools index ${SAMPLE}.protozoa.bam
samtools fastq ${SAMPLE}.virus.bam >${SAMPLE}.virus.fastq
samtools fastq ${SAMPLE}.fungi.bam >${SAMPLE}.fungi.fastq
samtools fastq ${SAMPLE}.protozoa.bam >${SAMPLE}.protozoa.fastq

$megahit --k-step 10 -r ${SAMPLE}.virus.fastq -o ${SAMPLE}.virus.megahit --out-prefix virus.megahit
$megahit --k-step 10 -r ${SAMPLE}.fungi.fastq -o ${SAMPLE}.fungi.megahit --out-prefix fungi.megahit
$megahit --k-step 10 -r ${SAMPLE}.protozoa.fastq -o ${SAMPLE}.protozoa.megahit --out-prefix protozoa.megahit
mv ${SAMPLE}.virus.megahit/virus.megahit.contigs.fa ${SAMPLE}.virus.megahit.contigs.fa
mv ${SAMPLE}.virus.megahit/fungi.megahit.contigs.fa ${SAMPLE}.fungi.megahit.contigs.fa
mv ${SAMPLE}.virus.megahit/protozoa.megahit.contigs.fa ${SAMPLE}.protozoa.megahit.contigs.fa



bwa index ${SAMPLE}.virus.megahit.contigs.fa
bwa mem  ${SAMPLE}.virus.megahit.contigs.fa ${SAMPLE}.virus.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.virus.bam


samtools depth ${SAMPLE}.megahit.contigs.virus.bam>${SAMPLE}.megahit.contigs.virus.cov
samtools view -H ${SAMPLE}.megahit.contigs.virus.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.virus.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.virus.uniq.cov


#fungi----
bwa index ${SAMPLE}.fungi.megahit.contigs.fa
bwa mem  ${SAMPLE}.fungi.megahit.contigs.fa ${SAMPLE}.fungi.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.fungi.bam
samtools depth ${SAMPLE}.megahit.contigs.fungi.bam>${SAMPLE}.megahit.contigs.fungi.cov
samtools view -H ${SAMPLE}.megahit.contigs.fungi.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.fungi.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.fungi.uniq.cov


#protozoa----
bwa index ${SAMPLE}.protozoa.megahit.contigs.fa
bwa mem  ${SAMPLE}.protozoa.megahit.contigs.fa ${SAMPLE}.protozoa.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.protozoa.bam
samtools depth ${SAMPLE}.megahit.contigs.protozoa.bam>${SAMPLE}.megahit.contigs.protozoa.cov
samtools view -H ${SAMPLE}.megahit.contigs.protozoa.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.protozoa.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.protozoa.uniq.cov



echo "--------3. Get rDNA coverage"

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




echo "--------4. Get MT coverage"

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


#This is only for Onco Paneles
#echo "--------5. TE elements"
#samtools bam2fq $BAM | bowtie2  -x ${DIR_CODE}/db.human/repeats/repbase.fa --end-to-end - | samtools view -F 4 -bh - | samtools sort - >${SAMPLE}.sort.repeat.bam


echo "--------5. Get TCR/BCR coverage"


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


echo "--------6. CNV"

readCounter  --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" $BAM > $SAMPLE.wig
PREFIX=$(basename $BAM | awk -F ".bam" '{print $1}')



Rscript ${DIR_CODE}/tools/ichorCNA/scripts/runIchorCNA.R --id $PREFIX \
--WIG $SAMPLE.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
--gcWig ${DIR_CODE}/db.human/ichorCNA/gc_hg19_1000kb.wig \
--mapWig ${DIR_CODE}/db.human/ichorCNA/map_hg19_1000kb.wig \
--centromere ${DIR_CODE}/db.human/ichorCNA/GRCh37.p13_centromere_UCSC-gapTable.txt \
--normalPanel ${DIR_CODE}/db.human/ichorCNA/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
--includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
--estimateNormal True --estimatePloidy True --estimateScPrevalence True \
--scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $OUTDIR/




echo "--------7. T and B cell repetoire profiling"

python /u/home/s/serghei/collab/code/imrep/imrep.py --bam --extendedOutput ${BAM} ${SAMPLE}.imrep.cdr3



rm -fr ${2}/*cleaned_input.fasta
rm -fr ${2}/*fastq
rm -fr ${2}/partial*


echo "--------8. calculate off target coverage"


while read line
do
chr=$(echo $line | awk '{print $1}')
x=$(echo $line | awk '{print $2}')
y=$(echo $line | awk '{print $3}')

n=0
n=$(samtools view -bh $BAM $chr:$x-$y | samtools  depth - | awk '{s+=$3} END {print s}')
echo $chr,$x,$y,$n >>${SAMPLE}.offtarget.cov
done<${DIR_CODE}/db.human/intergenic.regions/intergenic.regions.hg19.autosomes.bed

rm -fr ${2}/*bam
rm -fr ${2}/*bai


end=`date +%s`

runtime=$((end-start))

echo "It took "$runtime" to run SBT protocol on " $BAM





