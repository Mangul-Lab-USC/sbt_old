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

bam_mtDNA=${OUT}.mtDNA.sort.bam
header=${OUT}.header.txt
bam_mtDNA_unique=${OUT}.mtDNA.sort.unique.bam
cov=${OUT}.mtDNA.cov

rm -fr $OUT_DIR
mkdir $OUT_DIR


bowtie2  -x /PHShome/sv188/sbt/db.human/MT --end-to-end -1 $IN_FASTQ1 -2 $IN_FASTQ2 | samtools view -F 4 -bh - | samtools sort - >$bam_mtDNA




samtools index $bam_mtDNA

samtools view -H $bam_mtDNA >$header
samtools view -F 12  $bam_mtDNA | grep -v "XS:" | cat $header - | samtools view -b - > $bam_mtDNA_unique
samtools depth $bam_mtDNA_unique >$cov


rm -fr $bam_mtDNA
rm -fr ${bam_mtDNA}.bai
rm -fr $header
rm -fr $bam_mtDNA_unique



