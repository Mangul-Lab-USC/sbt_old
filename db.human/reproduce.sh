


for chr in {1..22}
do

awk -F "," '{if ($1=='$chr') print}' DNA.kmers.75.clean.cov.csv >rDNA.kmers.75.clean.cov.chr${chr}.csv
python infer.coverage.intervals.py rDNA.kmers.75.clean.cov.chr${chr}.csv ${chr} rDNA.kmers.75.clean.chr${chr}.bed
done




