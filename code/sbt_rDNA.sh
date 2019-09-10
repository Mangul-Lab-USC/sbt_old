#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('in_fastq1')
parser.add_argument('in_fastq2')
parser.add_argument('out_dir')
EOF

prefix=$(basename "$IN_FASTQ1" .bam)
OUT=$OUT_DIR"/"$prefix

mkdir $OUT_DIR
cd $OUT_DIR





module load samtools
module load bowtie2
module load bcftools

bam_rDNA=${OUT}.rDNA.sort.rDNA.bam
header=${OUT}.header.txt
bam_rDNA_unique=${OUT}.rDNA.sort.rDNA.unique.bam
cov=${OUT}.rDNA.cov
bcf=${OUT}.rDNA.bcf

rm -fr $OUT_DIR
mkdir $OUT_DIR


bowtie2  -x /PHShome/sv188/sbt/db.human/rDNA --end-to-end -1 $IN_FASTQ1 -2 $IN_FASTQ2 | samtools view -F 4 -bh - | samtools sort - >$bam_rDNA




samtools index $bam_rDNA

samtools view -H $bam_rDNA >$header
samtools view -F 12  $bam_rDNA | grep -v "XS:" | cat $header - | samtools view -b - > $bam_rDNA_unique
samtools depth $bam_rDNA_unique >$cov


rm -fr $bam_rDNA
rm -fr ${bam_rDNA}.bai
rm -fr $header
rm -fr $bam_rDNA_unique



