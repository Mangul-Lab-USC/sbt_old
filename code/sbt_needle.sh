#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('bam')
parser.add_argument('out_dir')
EOF

prefix=$(basename "$BAM" .bam)
OUT=$OUT_DIR"/"$prefix

mkdir $OUT_DIR





module load samtools
module load bowtie2
module load bcftools


mkdir ${OUT_DIR}_temp
cd ${OUT_DIR}_temp
pwd



/PHShome/sv188/needle/needle.sh -f ${BAM} ${OUT_DIR}_temp
mv ${OUT_DIR}_temp/*contigs.fa ${OUT_DIR}
mv ${OUT_DIR}_temp/*filtered.csv ${OUT_DIR}
mv ${OUT_DIR}_temp/*putative.bacteria.csv ${OUT_DIR}

rm -fr ${OUT_DIR}_temp


