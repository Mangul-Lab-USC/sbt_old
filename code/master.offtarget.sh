ls $PWD/*bam | awk -F ".bam" '{print $1}' >samples.txt

while read line
do

fbname=$(basename "$line" | cut -d. -f1)

echo "~/code/sbt/code/sbt_offtarget_cov.sh ${line}.bam $line">run.${fbname}.sh

#qsub -cwd -V -N off_coverage -l h_data=8G,highp,time=10:00:00 run.${line}.sh

done<samples.txt

~/code/miscellaneous.scripts/submit_QSUB_array.sh 1 10
