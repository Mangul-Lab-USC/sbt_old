Rscript prepare.intergenic.regions.R
awk '{print $2,$3,$4,$5}' intergenic.regions.hg19.txt | sed 's/chr//' | sed 's/\"//g' >intergenic.regions.hg19.bed
awk '{if ($1==1) print}' intergenic.regions.hg19.bed >intergenic.regions.hg19.chr1.bed

#autosomes
head -n 17555 intergenic.regions.hg19.bed >intergenic.regions.hg19.autosomes.bed



