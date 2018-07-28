library("GeneImp",lib="xxxxx")
arg = commandArgs(trailingOnly=T)
arg.vcf = arg[1]
arg.ref = arg[2]
imputevcf(vcfname=arg.vcf,ref.vcfname=arg.ref,maxjobs=8,klthresh=10,numfilterhaps=200,verbose=2)
