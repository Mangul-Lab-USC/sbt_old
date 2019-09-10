
ls *bam | awk -F ".bam" '{print $1}' >samples.txt

while read line
do

echo "~/code/sbt/code/sbt_offtarget_cov.sh ${line}.bam $line" >run.${line}.sh


#qsub -cwd -V -N imrep -l h_data=16G,highp,time=24:00:00 run.${line}.sh

done<samples.txt

~/code/miscellaneous.scripts/submit_QSUB_array.sh  16 24
