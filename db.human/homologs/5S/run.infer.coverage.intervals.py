for chr in {1..22}
do

sed 's/chr//' 5S.kmers.75bp.cov | awk '{print $1","$2","$3}' | awk -F "," '{if ($1=='$chr') print}' >5S.kmers.75.clean.cov.chr${chr}.csv
python infer.coverage.intervals.py 5S.kmers.75.clean.cov.chr${chr}.csv ${chr} 5S.kmers.75.clean.chr${chr}.bed
done
