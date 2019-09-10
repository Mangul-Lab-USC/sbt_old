ls *bam | awk -F ".bam" '{print $1}' >samples.txt

while read line
do


echo "samtools index ${line}.bam">run.${line}.sh

#qsub -cwd -V -N off_coverage -l h_data=8G,highp,time=10:00:00 run.${line}.sh

done<samples.txt

~/code/miscellaneous.scripts/submit_QSUB_array.sh 4 1
