

. /u/local/Modules/default/init/modules.sh
module load bwa
samtools=~/collab/code/rop.old/tools/samtools


bwa mem ~/project/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa rDNA.fasta | $samtools view -S -b -F 4 - | $samtools sort - >rDNA.bam

$samtools index rDNA.bam

$samtools depth rDNA.bam >rDNA.cov

