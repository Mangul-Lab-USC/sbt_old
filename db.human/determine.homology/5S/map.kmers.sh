. /u/local/Modules/default/init/modules.sh
module load bwa
samtools=~/collab/code/rop.old/tools/samtools


bwa mem /u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa 5S.kmers.75bp.fasta | $samtools view -S -b -F 4 - | $samtools sort - >5S.kmers.75bp.bam

$samtools index 5S.kmers.75bp.bam

$samtools depth 5S.kmers.75bp.bam >5S.kmers.75bp.cov
