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

start=`date +%s`

echo  "Start. "$start

. /u/local/Modules/default/init/modules.sh
module load bowtie2

#samtools=/PHShome/sv188/code/rop/tools/samtools-1.3/samtools
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
bam=$1
prefix=$(basename $bam | awk -F ".bam" '{print $1}')

DIR=$2
mkdir $DIR

SAMPLE=$DIR"/"${prefix}



DIR_CODE=$(dirname $0)
megahit=/u/home/n/nlapier2/lib/megahit_v1.1.3_LINUX_CPUONLY_x86_64-bin/megahit
metaphlan=/u/home/s/serghei//collab/code/rop/tools/metaphlan2/metaphlan2.py
imrep=/u/home/s/serghei//collab/code/imrep/imrep.py


#0 cov of intergenic
rm -fr ${SAMPLE}.offtarget.cov





#1. Extract unmapped reads
echo "---------------------------------------------------------------------------"
echo "Extract unmapped reads from " $1
echo "$samtools view -f 0x4 -bh $bam | $samtools bam2fq - >${SAMPLE}.unmapped.fastq"
$samtools view -f 0x4 -bh $bam | $samtools bam2fq - >${SAMPLE}.unmapped.fastq






#repeat
#rm -fr ${SAMPLE}.repeat.fastq
#while read line;do chr=$(echo $line | awk -F "," '{print $1}');x=$(echo $line | awk -F "," '{print $2}');y=$(echo $line | awk -F "," '{print $3}');$samtools view -bh $bam $chr:$x-$y | $samtools fastq - >>${SAMPLE}.repeat.fastq;done<${DIR_CODE}/db.human/repeats/repbase.kmers.75.clean.filtered.bed 
#cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.repeat.fastq >${SAMPLE}.cat.repeat.fastq 
#rm -fr ${SAMPLE}.repeat.fastq
#bowtie2  -x ${DIR_CODE}/db.human/repeats/repbase.fa --end-to-end ${SAMPLE}.cat.repeat.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.repeat.bam
#$samtools index ${SAMPLE}.sort.repeat.bam
#$samtools fastq ${SAMPLE}.sort.repeat.bam>${SAMPLE}.sort.repeat.fastq
#${DIR_CODE}/megahit/megahit --k-step 10 -r ${SAMPLE}.sort.repeat.fastq -o ${SAMPLE}.repeat.megahit --out-prefix repeat.megahit
#$samtools mpileup -f  ${DIR_CODE}/db.human/repeats/repbase.fa ${SAMPLE}.sort.repeat.bam >${SAMPLE}.repeat.pileup
#$samtools mpileup -uf ${DIR_CODE}/db.human/repeats/repbase.fa ${SAMPLE}.sort.repeat.bam | bcftools  call -mv -Oz >${SAMPLE}.repeat.bcf


#map all reads onto TE----
#$samtools bam2fq $bam >${SAMPLE}.fastq


#bowtie2  -x ${DIR_CODE}/db.human/repeats/repbase.fa --end-to-end ${SAMPLE}.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.repeat.bam
#$samtools mpileup -f  ${DIR_CODE}/db.human/repeats/repbase.fa ${SAMPLE}.sort.repeat.bam >${SAMPLE}.repeat.pileup
#$samtools mpileup -uf ${DIR_CODE}/db.human/repeats/repbase.fa ${SAMPLE}.sort.repeat.bam | bcftools  call -mv -Oz >${SAMPLE}.repeat.bcf

#rm -fr ${SAMPLE}.fastq

#megahit
$samtools view -bh $bam NC_007605 | $samtools fastq - > ${SAMPLE}.NC_007605.fastq
cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.NC_007605.fastq>${SAMPLE}.cat.unmapped.fastq

$megahit --k-step 10 -r ${SAMPLE}.cat.unmapped.fastq -o ${SAMPLE}.unmapped.megahit --out-prefix ${prefix}.unmapped.megahit


mv ${SAMPLE}.unmapped.megahit/${prefix}.unmapped.megahit.contigs.fa ${DIR}
gzip ${prefix}.unmapped.megahit.contigs.fa
rm -fr ${SAMPLE}.unmapped.megahit



module load bwa 
bwa mem -a ${DIR_CODE}/db.human/viral.vipr/NONFLU_All.fastq ${SAMPLE}.cat.unmapped.fastq | $samtools view -S -b -F 4 - | $samtools sort - >${SAMPLE}.virus.bam
bwa mem -a ${DIR_CODE}/db.human/fungi/fungi.ncbi.february.3.2018.fasta ${SAMPLE}.unmapped.fastq | $samtools view -S -b -F 4 - |  $samtools sort - >${SAMPLE}.fungi.bam
bwa mem -a ${DIR_CODE}/db.human/protozoa/protozoa.ncbi.february.3.2018.fasta ${SAMPLE}.unmapped.fastq | $samtools view -S -b -F 4 - | $samtools sort - >${SAMPLE}.protozoa.bam

$samtools index ${SAMPLE}.virus.bam
$samtools index ${SAMPLE}.fungi.bam 
$samtools index ${SAMPLE}.protozoa.bam

$samtools mpileup -f ${DIR_CODE}/db.human/viral.vipr/NONFLU_All.fastq ${SAMPLE}.virus.bam | gzip - >${SAMPLE}.virus.pileup.gz

$samtools mpileup -f ${DIR_CODE}/db.human/fungi/fungi.ncbi.february.3.2018.fasta ${SAMPLE}.fungi.bam | gzip - >${SAMPLE}.fungi.pileup.gz

$samtools mpileup -f ${DIR_CODE}/db.human/protozoa/protozoa.ncbi.february.3.2018.fasta ${SAMPLE}.protozoa.bam | gzip - >${SAMPLE}.protozoa.pileup.gz 



$samtools fastq ${SAMPLE}.virus.bam >${SAMPLE}.virus.fastq
$samtools fastq ${SAMPLE}.fungi.bam >${SAMPLE}.fungi.fastq
$samtools fastq ${SAMPLE}.protozoa.bam >${SAMPLE}.protozoa.fastq

$megahit --k-step 10 -r ${SAMPLE}.virus.fastq -o ${SAMPLE}.virus.megahit --out-prefix ${prefix}.virus.megahit
mv ${SAMPLE}.virus.megahit/${prefix}.virus.megahit.contigs.fa ${DIR}
gzip ${prefix}.virus.megahit.contigs.fa
rm -fr ${SAMPLE}.virus.megahit


$megahit --k-step 10 -r ${SAMPLE}.fungi.fastq -o ${SAMPLE}.fungi.megahit --out-prefix ${prefix}.fungi.megahit
mv ${SAMPLE}.fungi.megahit/${prefix}.fungi.megahit.contigs.fa ${DIR}
gzip ${prefix}.fungi.megahit.contigs.fa
rm -fr ${SAMPLE}.fungi.megahit


$megahit --k-step 10 -r ${SAMPLE}.protozoa.fastq -o ${SAMPLE}.protozoa.megahit --out-prefix ${prefix}.protozoa.megahit
mv ${SAMPLE}.protozoa.megahit/${prefix}.protozoa.megahit.contigs.fa ${DIR}
gzip ${prefix}.protozoa.megahit.contigs.fa
rm -fr ${SAMPLE}.protozoa.megahit


#2. Extract reads from GL000220.1
echo "---------------------------------------------------------------------------"
echo "Extract reads from GL000220.1 and save to "${SAMPLE}.GL000220.fastq
$samtools view -bh ${bam} GL000220.1 | $samtools bam2fq - >${SAMPLE}.GL000220.fastq

rm -fr ${SAMPLE}.rDNA.mapped.fastq
while read line;do chr=$(echo $line | awk -F "," '{print $1}');x=$(echo $line | awk -F "," '{print $2}');y=$(echo $line | awk -F "," '{print $3}');$samtools view -bh $bam $chr:$x-$y | $samtools fastq - >>${SAMPLE}.rDNA.mapped.fastq;done<${DIR_CODE}/db.human/rDNA.kmers.75.clean.filtered.bed




# 4. merge unmapped + GL000220.1 + 5S and map to 45S and 5S
echo "---------------------------------------------------------------------------"
echo "Merge unmapped + GL000220.1 + 5S and map to 45S and 5S " ${SAMPLE}.cat.rDNA.fastq
module load bcftools

cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.GL000220.fastq ${SAMPLE}.rDNA.mapped.fastq >${SAMPLE}.cat.rDNA.fastq


bowtie2  -x ${DIR_CODE}/db.human/rDNA --end-to-end ${SAMPLE}.cat.rDNA.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.rDNA.bam

$samtools index ${SAMPLE}.sort.rDNA.bam


$samtools mpileup -f  ${DIR_CODE}/db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.bam | gzip - >${SAMPLE}.rDNA.pileup.gzip
$samtools mpileup -uf ${DIR_CODE}/db.human/rDNA.fasta ${SAMPLE}.sort.rDNA.bam | bcftools  call -mv -Oz >${SAMPLE}.rDNA.bcf

#$samtools fastq ${SAMPLE}.sort.rDNA.bam>${SAMPLE}.sort.rDNA.fastq
#${DIR_CODE}/megahit/megahit --k-step 10 -r ${SAMPLE}.sort.rDNA.fastq -o ${SAMPLE}.rDNA.megahit --out-prefix rDNA.megahit

rm -fr ${SAMPLE}.sort.rDNA.bam
rm -fr ${SAMPLE}.cat.rDNA.fastq ${SAMPLE}.GL000220.fastq ${SAMPLE}.rDNA.mapped.fastq ${SAMPLE}.sort.rDNA.fastq 


# --------------------------------------------------------
# 5.merge MT and unmapped and map to MT
echo "---------------------------------------------------------------------------"
echo "merge MT and unmapped and map to MT"

$samtools view -bh $bam MT | $samtools bam2fq - >${SAMPLE}.MT.fastq
cat ${SAMPLE}.MT.fastq ${SAMPLE}.unmapped.fastq >${SAMPLE}.cat.MT.fastq

#${DIR_CODE}/megahit/megahit --k-step 10 -r ${SAMPLE}.cat.MT.fastq -o ${SAMPLE}.MT.megahit --out-prefix MT.megahit
bowtie2  -x ${DIR_CODE}/db.human/MT --end-to-end ${SAMPLE}.cat.MT.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.MT.bam

$samtools index ${SAMPLE}.sort.MT.bam

$samtools fastq ${SAMPLE}.sort.MT.bam >${SAMPLE}.sort.MT.fastq

#${DIR_CODE}/megahit/megahit --k-step 10 -r ${SAMPLE}.sort.MT.fastq -o ${SAMPLE}.MT.megahit --out-prefix MT.megahit


$samtools mpileup -f ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.sort.MT.bam | gzip - >${SAMPLE}.MT.pileup.gzip
$samtools mpileup -uf ${DIR_CODE}/db.human/MT.fasta ${SAMPLE}.sort.MT.bam | bcftools  call -mv -Oz >${SAMPLE}.MT.bcf


rm -fr ${SAMPLE}.sort.MT.bam
rm -fr ${SAMPLE}.cat.MT.fastq ${SAMPLE}.sort.MT.fastq
rm -fr ${SAMPLE}.MT.fastq


#TODO - coverage of off target regions - background coverage 




# --------------------------------------------------------
#6 extrcat TCR/BCR reads run imrep
echo "---------------------------------------------------------------------------"
echo "extract immune reads and map to VDJ ref"
#IGH
$samtools view -bh $bam 14:106032614-107288051 |  $samtools bam2fq - >${SAMPLE}.ireceptor.fastq
#IGK
$samtools view -bh $bam 2:89156874-89630436 |  $samtools bam2fq - >>${SAMPLE}.ireceptor.fastq
#IGL
$samtools view -bh $bam 22:22380474-23265085 |  $samtools bam2fq - >>${SAMPLE}.ireceptor.fastq

#TCRA
$samtools view -bh $bam 14:22090057-23021075 |  $samtools bam2fq - >>${SAMPLE}.ireceptor.fastq

#TCRB
$samtools view -bh $bam 7:141998851-141998851 |  $samtools bam2fq - >>${SAMPLE}.ireceptor.fastq

#TCRG
$samtools view -bh $bam 7:38279625-38407656 |  $samtools bam2fq - >>${SAMPLE}.ireceptor.fastq

cat ${SAMPLE}.ireceptor.fastq ${SAMPLE}.unmapped.fastq >${SAMPLE}.cat.ireceptor.fastq

bowtie2  -x ${DIR_CODE}/db.human/IMGT/receptor.genes ${SAMPLE}.cat.ireceptor.fastq | $samtools view -F 4 -bh - | $samtools sort - >${SAMPLE}.sort.ireceptor.bam

$samtools index ${SAMPLE}.sort.ireceptor.bam

#$samtools fastq ${SAMPLE}.sort.ireceptor.bam >${SAMPLE}.sort.VDJ.fastq
#${DIR_CODE}/megahit/megahit --k-min 51 --k-step 10 -r ${SAMPLE}.sort.VDJ.fastq -o ${SAMPLE}.VDJ.megahit --out-prefix VDJ.megahit




$samtools mpileup -f ${DIR_CODE}/db.human/IMGT/receptor.genes.fasta ${SAMPLE}.sort.ireceptor.bam | gzip - >${SAMPLE}.ireceptor.pileup.gzip
$samtools mpileup -uf ${DIR_CODE}/db.human/IMGT/receptor.genes.fasta ${SAMPLE}.sort.ireceptor.bam | bcftools  call -mv -Oz >${SAMPLE}.ireceptor.bcf


rm -fr ${SAMPLE}.sort.ireceptor.bam
rm -fr ${SAMPLE}.ireceptor.fastq
rm -fr ${SAMPLE}.cat.ireceptor.fastq


#7 imrep
echo "---------------------------------------------------------------------------"
# 8 map to microbiome - sort of micop here! do better fungi reference TODO!
#metaphlan

# assemle contigs and map to viruses
echo "---------------------------------------------------------------------------"








module load bwa




python $metaphlan ${SAMPLE}.unmapped.fastq --nproc 8 --input_type fastq --bowtie2out ${SAMPLE}.metaphlan2.bowtie2.txt >${SAMPLE}.metaphlan2.txt




  
python $imrep --bam --extendedOutput ${bam} ${SAMPLE}.imrep.cdr3



rm -fr ${SAMPLE}.unmapped.fastq
rm -fr ${2}/*cleaned_input.fasta
rm -fr ${2}/*fastq
rm -fr ${2}/partial*


echo "0. Get off target coverage"

while read line
do
chr=$(echo $line | awk '{print $1}')
x=$(echo $line | awk '{print $2}')
y=$(echo $line | awk '{print $3}')

n=0
n=$($samtools view -bh $bam $chr:$x-$y | $samtools  depth - | awk '{s+=$3} END {print s}')
echo $chr,$x,$y,$n >>${SAMPLE}.offtarget.cov
done<${DIR_CODE}/db.human/intergenic.regions/intergenic.regions.hg19.autosomes.bed

rm -fr ${2}/*bam 
rm -fr ${2}/*bai


end=`date +%s`

runtime=$((end-start))

echo "It took "$runtime" to run SBT protocol on " $bam
