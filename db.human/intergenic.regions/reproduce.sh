Rscript prepare.intergenic.regions.R
awk '{print $2,$3,$4,$5}' intergenic.regions.hg19.txt | sed 's/chr//' | sed 's/\"//g' >intergenic.regions.hg19.bed

