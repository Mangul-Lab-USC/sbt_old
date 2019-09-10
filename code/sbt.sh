#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('bam')
parser.add_argument('out_dir')
EOF

rm -fr $OUT_DIR
mkdir $OUT_DIR

cd $OUT_DIR

PREFIX=$(basename "$BAM" .bam)

FASTQ_1=${OUT_DIR}/${PREFIX}_1.fastq
FASTQ_2=${OUT_DIR}/${PREFIX}_2.fastq


echo "module load samtools">master_${PREFIX}.sh

echo "samtools fastq $BAM -1 $FASTQ_1 -2 $FASTQ_2">>master_${PREFIX}.sh

echo "/PHShome/sv188/sbt/sbt_rDNA.sh $FASTQ_1 $FASTQ_2 ${OUT_DIR}/${PREFIX}_rDNA/">>master_${PREFIX}.sh
echo "/PHShome/sv188/sbt/sbt_mtDNA.sh $FASTQ_1 $FASTQ_2 ${OUT_DIR}/${PREFIX}_mtDNA/">>master_${PREFIX}.sh





echo "/PHShome/sv188/sbt/sbt_TE.sh $FASTQ_1 $FASTQ_2 ${OUT_DIR}/${PREFIX}_TE/">>master_${PREFIX}.sh

echo "/PHShome/sv188/sbt/sbt_needle.sh ${BAM} ${OUT_DIR}/${PREFIX}_needle/">>master_${PREFIX}.sh
echo "/PHShome/sv188/sbt/sbt_offtarget_cov.sh ${BAM} ${OUT_DIR}/${PREFIX}_offcov/">>master_${PREFIX}.sh
echo "/PHShome/sv188/sbt/sbt_imrep.sh ${BAM} ${OUT_DIR}/${PREFIX}_imrep/">>master_${PREFIX}.sh

chmod 755 master_${PREFIX}.sh 
./master_${PREFIX}.sh 

rm -fr $FASTQ_1 $FASTQ_2

