. /u/local/Modules/default/init/modules.sh
module load bwa
samtools=~/collab/code/rop.old/tools/samtools


bwa mem /u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa HSU13369.kmers.75bp.fasta | $samtools view -S -b -F 4 - | $samtools sort - >HSU13369.kmers.75bp.bam

$samtools index HSU13369.kmers.75bp.bam

$samtools depth HSU13369.kmers.75bp.bam >HSU13369.kmers.75bp.cov
