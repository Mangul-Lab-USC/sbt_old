# seeing.beyond.target
Seeing Beyond the Target (SBT) - computational protocol to process off target reads in WXS experiments

-- Error

```
~/project/code/seeing.beyond.target/sbt.sh -op LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.bam test
Note: Assumptions about bam file : 1. mapped to hg19 2. mitochondria tag is MT
infile: 
outfile: test
Add MiniConda to PATH if it's available
Start SBT analysis ... 1554019501
------------- ... Extract unmapped reads from ... /u/project/zarlab/serghei/op_test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.bam -------------
Number of non human reference in the BAM
1 /u/project/zarlab/serghei/op_test/test/non.human.references.txt
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 2221423 reads
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 3749 reads
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 15535 reads
62140 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 432 reads
63868 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 16523 reads
129960 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 43377 reads
303468 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 0 reads
303468 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 16181 reads
368192 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.ireceptor.fastq
Parse bam file with mapped and unmapped reads
--->chr 14
('Number of reads extacted from ', 'IGH', 'locus : ', 5049)
('Number of reads extacted from ', 'IGK', 'locus : ', 125)
('Number of reads extacted from ', 'IGL', 'locus : ', 5450)
('Number of reads extacted from ', 'TRA', 'locus : ', 15517)
('Number of reads extacted from ', 'TRB', 'locus : ', 0)
('Number of reads extacted from ', 'TRG', 'locus : ', 4801)
('Number of unmapped reads extracted', 2221423)
Starting ImReP-0.6 (developped by Igor Mandric, Serghei Mangul)


 _                          
(_)_ __ ___  _ __ ___ _ __  
| | '_ ` _ \| '__/ _ \ '_ \ 
| | | | | | | | |  __/ |_) |
|_|_| |_| |_|_|  \___| .__/ 
                     |_|    




hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
hello
113518 partial-V CDR3 found
6838 partial-J CDR3 found
90 full CDR3 found:
	- 1 of type IGH
	- 0 of type IGK
	- 0 of type IGL
	- 88 of type TRA
	- 0 of type TRB
	- 0 of type TRD
	- 1 of type TRG
Done. Bye-bye
NC_007605
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 0 reads
------------- (3) metaphlan2 -------------
wc: /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.metaphlan2.txt: No such file or directory
------------- (3) Use needle to detect viruses, fungia,and protozoa -------------
required infile: 
required outfile: /u/project/zarlab/serghei/op_test/test/needle.out
Start needle analysis ... 1554020537
Extract unmapped reads from  /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq
FASTQ is provided
/u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/needle.sh: line 105: ------>Number of unmapped reads: command not found
8885692 /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 99010 sequences (10000010 bp)...
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 15.220 CPU sec, 15.108 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 16.118 CPU sec, 15.963 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 15.341 CPU sec, 15.165 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 16.396 CPU sec, 16.245 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.608 CPU sec, 19.491 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.952 CPU sec, 19.797 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.453 CPU sec, 19.286 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.252 CPU sec, 19.099 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.530 CPU sec, 19.358 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.875 CPU sec, 19.738 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.534 CPU sec, 19.360 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.289 CPU sec, 19.132 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.509 CPU sec, 19.340 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 20.389 CPU sec, 20.240 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 20.211 CPU sec, 20.043 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.583 CPU sec, 19.434 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.498 CPU sec, 19.334 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.450 CPU sec, 19.294 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.346 CPU sec, 19.178 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.986 CPU sec, 19.853 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.807 CPU sec, 19.664 real sec
[M::process] read 43203 sequences (4363503 bp)...
[M::mem_process_seqs] Processed 99010 reads in 19.979 CPU sec, 19.865 real sec
[M::mem_process_seqs] Processed 43203 reads in 8.190 CPU sec, 8.099 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/viral.vipr/NONFLU_All.fastq /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq
[main] Real time: 434.429 sec; CPU: 433.368 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 99010 sequences (10000010 bp)...
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 10.944 CPU sec, 10.836 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 11.279 CPU sec, 11.135 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 10.990 CPU sec, 10.836 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 11.689 CPU sec, 11.556 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.737 CPU sec, 13.588 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.993 CPU sec, 13.855 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.973 CPU sec, 13.823 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.880 CPU sec, 13.741 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.991 CPU sec, 13.842 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.030 CPU sec, 13.885 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.950 CPU sec, 13.798 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.929 CPU sec, 13.789 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.061 CPU sec, 13.916 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.343 CPU sec, 14.203 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.378 CPU sec, 14.228 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.089 CPU sec, 13.961 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 13.848 CPU sec, 13.698 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.109 CPU sec, 13.977 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.018 CPU sec, 13.878 real sec
[M::process] read 99010 sequences (10000010 bp)...



[M::mem_process_seqs] Processed 99010 reads in 14.502 CPU sec, 14.369 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.930 CPU sec, 14.777 real sec
[M::process] read 43203 sequences (4363503 bp)...
[M::mem_process_seqs] Processed 99010 reads in 14.090 CPU sec, 13.989 real sec
[M::mem_process_seqs] Processed 43203 reads in 5.706 CPU sec, 5.627 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/fungi/fungi.ncbi.february.3.2018.fasta /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq
[main] Real time: 304.025 sec; CPU: 306.240 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 99010 sequences (10000010 bp)...
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 8.617 CPU sec, 8.507 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 8.909 CPU sec, 8.764 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 8.684 CPU sec, 8.532 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 8.741 CPU sec, 8.609 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.417 CPU sec, 9.267 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.635 CPU sec, 9.501 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.364 CPU sec, 9.207 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.593 CPU sec, 9.450 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.378 CPU sec, 9.231 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.647 CPU sec, 9.512 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.335 CPU sec, 9.185 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.430 CPU sec, 9.298 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.445 CPU sec, 9.294 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.662 CPU sec, 9.523 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.703 CPU sec, 9.563 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.493 CPU sec, 9.352 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.371 CPU sec, 9.218 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.472 CPU sec, 9.335 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.379 CPU sec, 9.227 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.338 CPU sec, 9.203 real sec
[M::process] read 99010 sequences (10000010 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.441 CPU sec, 9.295 real sec
[M::process] read 43203 sequences (4363503 bp)...
[M::mem_process_seqs] Processed 99010 reads in 9.302 CPU sec, 9.197 real sec
[M::mem_process_seqs] Processed 43203 reads in 3.878 CPU sec, 3.802 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/protozoa/protozoa.ncbi.february.3.2018.fasta /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq
[main] Real time: 207.041 sec; CPU: 209.677 sec
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 389670 reads
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 62528 reads
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 44955 reads
23.0Gb memory in total.
Using: 21.234Gb.
MEGAHIT v1.1.3
--- [Sun Mar 31 01:38:23 2019] Start assembly. Number of CPU threads 12 ---
--- [Sun Mar 31 01:38:23 2019] Available memory: 25333133312, used: 22799819980
--- [Sun Mar 31 01:38:23 2019] Converting reads to binaries ---
    [read_lib_functions-inl.h  : 209]     Lib 0 (/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.fastq): se, 389421 reads, 101 max length
    [utils.h                   : 126]     Real: 0.4576	user: 0.2850	sys: 0.0650	maxrss: 24692
--- [Sun Mar 31 01:38:23 2019] k list: 21,31,41,51,61,71,81,91,101,111,121,131,141 ---
--- [Sun Mar 31 01:38:23 2019] Extracting solid (k+1)-mers for k = 21 ---
--- [Sun Mar 31 01:38:25 2019] Building graph for k = 21 ---
--- [Sun Mar 31 01:38:27 2019] Assembling contigs from SdBG for k = 21 ---
--- [Sun Mar 31 01:38:31 2019] Local assembling k = 21 ---
--- [Sun Mar 31 01:38:31 2019] Extracting iterative edges from k = 21 to 31 ---
--- [Sun Mar 31 01:38:32 2019] Building graph for k = 31 ---
--- [Sun Mar 31 01:38:33 2019] Assembling contigs from SdBG for k = 31 ---
--- [Sun Mar 31 01:38:34 2019] Local assembling k = 31 ---
--- [Sun Mar 31 01:38:35 2019] Extracting iterative edges from k = 31 to 41 ---
--- [Sun Mar 31 01:38:35 2019] Building graph for k = 41 ---
--- [Sun Mar 31 01:38:36 2019] Assembling contigs from SdBG for k = 41 ---
--- [Sun Mar 31 01:38:36 2019] Local assembling k = 41 ---
--- [Sun Mar 31 01:38:36 2019] Extracting iterative edges from k = 41 to 51 ---
--- [Sun Mar 31 01:38:37 2019] Building graph for k = 51 ---
--- [Sun Mar 31 01:38:37 2019] Assembling contigs from SdBG for k = 51 ---
--- [Sun Mar 31 01:38:37 2019] Local assembling k = 51 ---
--- [Sun Mar 31 01:38:37 2019] Extracting iterative edges from k = 51 to 61 ---
--- [Sun Mar 31 01:38:38 2019] Building graph for k = 61 ---
--- [Sun Mar 31 01:38:38 2019] Assembling contigs from SdBG for k = 61 ---
--- [Sun Mar 31 01:38:38 2019] Local assembling k = 61 ---
--- [Sun Mar 31 01:38:38 2019] Extracting iterative edges from k = 61 to 71 ---
--- [Sun Mar 31 01:38:39 2019] Building graph for k = 71 ---
--- [Sun Mar 31 01:38:39 2019] Assembling contigs from SdBG for k = 71 ---
--- [Sun Mar 31 01:38:39 2019] Local assembling k = 71 ---
--- [Sun Mar 31 01:38:40 2019] Extracting iterative edges from k = 71 to 81 ---
--- [Sun Mar 31 01:38:40 2019] Building graph for k = 81 ---
--- [Sun Mar 31 01:38:40 2019] Assembling contigs from SdBG for k = 81 ---
--- [Sun Mar 31 01:38:40 2019] Local assembling k = 81 ---
--- [Sun Mar 31 01:38:41 2019] Extracting iterative edges from k = 81 to 91 ---
--- [Sun Mar 31 01:38:41 2019] Merging to output final contigs ---
--- [STAT] 0 contigs, total 0 bp, min 0 bp, max 0 bp, avg 0 bp, N50 0 bp
--- [Sun Mar 31 01:38:41 2019] ALL DONE. Time elapsed: 18.465353 seconds ---
23.0Gb memory in total.
Using: 21.234Gb.
MEGAHIT v1.1.3
--- [Sun Mar 31 01:38:41 2019] Start assembly. Number of CPU threads 12 ---
--- [Sun Mar 31 01:38:41 2019] Available memory: 25333133312, used: 22799819980
--- [Sun Mar 31 01:38:41 2019] Converting reads to binaries ---
    [read_lib_functions-inl.h  : 209]     Lib 0 (/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.fastq): se, 61837 reads, 101 max length
    [utils.h                   : 126]     Real: 0.0946	user: 0.0490	sys: 0.0090	maxrss: 6348
--- [Sun Mar 31 01:38:41 2019] k list: 21,31,41,51,61,71,81,91,101,111,121,131,141 ---
--- [Sun Mar 31 01:38:41 2019] Extracting solid (k+1)-mers for k = 21 ---
--- [Sun Mar 31 01:38:42 2019] Building graph for k = 21 ---
--- [Sun Mar 31 01:38:43 2019] Assembling contigs from SdBG for k = 21 ---
--- [Sun Mar 31 01:38:44 2019] Local assembling k = 21 ---
--- [Sun Mar 31 01:38:44 2019] Extracting iterative edges from k = 21 to 31 ---
--- [Sun Mar 31 01:38:44 2019] Building graph for k = 31 ---
--- [Sun Mar 31 01:38:45 2019] Assembling contigs from SdBG for k = 31 ---
--- [Sun Mar 31 01:38:45 2019] Local assembling k = 31 ---
--- [Sun Mar 31 01:38:46 2019] Extracting iterative edges from k = 31 to 41 ---
--- [Sun Mar 31 01:38:46 2019] Building graph for k = 41 ---
--- [Sun Mar 31 01:38:46 2019] Assembling contigs from SdBG for k = 41 ---
--- [Sun Mar 31 01:38:46 2019] Local assembling k = 41 ---
--- [Sun Mar 31 01:38:46 2019] Extracting iterative edges from k = 41 to 51 ---
--- [Sun Mar 31 01:38:47 2019] Building graph for k = 51 ---
--- [Sun Mar 31 01:38:47 2019] Assembling contigs from SdBG for k = 51 ---
--- [Sun Mar 31 01:38:47 2019] Local assembling k = 51 ---
--- [Sun Mar 31 01:38:47 2019] Extracting iterative edges from k = 51 to 61 ---
--- [Sun Mar 31 01:38:47 2019] Building graph for k = 61 ---
--- [Sun Mar 31 01:38:47 2019] Assembling contigs from SdBG for k = 61 ---
--- [Sun Mar 31 01:38:47 2019] Local assembling k = 61 ---
--- [Sun Mar 31 01:38:47 2019] Extracting iterative edges from k = 61 to 71 ---
--- [Sun Mar 31 01:38:47 2019] Merging to output final contigs ---
--- [STAT] 2 contigs, total 435 bp, min 203 bp, max 232 bp, avg 218 bp, N50 232 bp
--- [Sun Mar 31 01:38:48 2019] ALL DONE. Time elapsed: 6.468047 seconds ---
23.0Gb memory in total.
Using: 21.234Gb.
MEGAHIT v1.1.3
--- [Sun Mar 31 01:38:48 2019] Start assembly. Number of CPU threads 12 ---
--- [Sun Mar 31 01:38:48 2019] Available memory: 25333133312, used: 22799819980
--- [Sun Mar 31 01:38:48 2019] Converting reads to binaries ---
    [read_lib_functions-inl.h  : 209]     Lib 0 (/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.fastq): se, 43598 reads, 101 max length
    [utils.h                   : 126]     Real: 0.0758	user: 0.0340	sys: 0.0090	maxrss: 6348
--- [Sun Mar 31 01:38:48 2019] k list: 21,31,41,51,61,71,81,91,101,111,121,131,141 ---
--- [Sun Mar 31 01:38:48 2019] Extracting solid (k+1)-mers for k = 21 ---
--- [Sun Mar 31 01:38:48 2019] Building graph for k = 21 ---
--- [Sun Mar 31 01:38:49 2019] Assembling contigs from SdBG for k = 21 ---
--- [Sun Mar 31 01:38:50 2019] Local assembling k = 21 ---
--- [Sun Mar 31 01:38:50 2019] Extracting iterative edges from k = 21 to 31 ---
--- [Sun Mar 31 01:38:50 2019] Building graph for k = 31 ---
--- [Sun Mar 31 01:38:51 2019] Assembling contigs from SdBG for k = 31 ---
--- [Sun Mar 31 01:38:51 2019] Local assembling k = 31 ---
--- [Sun Mar 31 01:38:51 2019] Extracting iterative edges from k = 31 to 41 ---
--- [Sun Mar 31 01:38:52 2019] Building graph for k = 41 ---
--- [Sun Mar 31 01:38:52 2019] Assembling contigs from SdBG for k = 41 ---
--- [Sun Mar 31 01:38:52 2019] Local assembling k = 41 ---
--- [Sun Mar 31 01:38:52 2019] Extracting iterative edges from k = 41 to 51 ---
--- [Sun Mar 31 01:38:52 2019] Building graph for k = 51 ---
--- [Sun Mar 31 01:38:53 2019] Assembling contigs from SdBG for k = 51 ---
--- [Sun Mar 31 01:38:53 2019] Local assembling k = 51 ---
--- [Sun Mar 31 01:38:53 2019] Extracting iterative edges from k = 51 to 61 ---
--- [Sun Mar 31 01:38:53 2019] Building graph for k = 61 ---
--- [Sun Mar 31 01:38:53 2019] Assembling contigs from SdBG for k = 61 ---
--- [Sun Mar 31 01:38:53 2019] Local assembling k = 61 ---
--- [Sun Mar 31 01:38:53 2019] Extracting iterative edges from k = 61 to 71 ---
--- [Sun Mar 31 01:38:53 2019] Merging to output final contigs ---
--- [STAT] 1 contigs, total 232 bp, min 232 bp, max 232 bp, avg 232 bp, N50 232 bp
--- [Sun Mar 31 01:38:53 2019] ALL DONE. Time elapsed: 5.656186 seconds ---
/u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/needle.sh: line 132: -------->Number of assembled contigs: command not found
0 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.fa
4 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.fa
2 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.fa
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.fa
[main] Real time: 0.035 sec; CPU: 0.005 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 99010 sequences (10000010 bp)...
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.fa
[main] Real time: 0.037 sec; CPU: 0.005 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 61837 sequences (6245537 bp)...
[M::mem_process_seqs] Processed 61837 reads in 1.503 CPU sec, 1.505 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.fa /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.fastq
[main] Real time: 1.747 sec; CPU: 1.594 sec
[bwa_index] Pack FASTA... 0.00 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.00 seconds elapse.
[bwa_index] Update BWT... 0.00 sec
[bwa_index] Pack forward-only FASTA... 0.00 sec
[bwa_index] Construct SA from BWT and Occ... 0.00 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.fa
[main] Real time: 0.034 sec; CPU: 0.005 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 43598 sequences (4403398 bp)...
[M::mem_process_seqs] Processed 43598 reads in 0.914 CPU sec, 0.914 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.fa /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.fastq
[main] Real time: 1.077 sec; CPU: 0.978 sec
-----------------------------------------------------
Map assembled contigs onto the microbial references
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/viral.vipr/NONFLU_All.fastq /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.fa
[main] Real time: 6.429 sec; CPU: 6.279 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 2 sequences (435 bp)...
[M::mem_process_seqs] Processed 2 reads in 0.008 CPU sec, 0.007 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/fungi/fungi.ncbi.february.3.2018.fasta /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.fa
[main] Real time: 0.661 sec; CPU: 0.663 sec
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 1 sequences (232 bp)...
[M::mem_process_seqs] Processed 1 reads in 0.002 CPU sec, 0.001 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -a /u/project/zarlab/serghei/code/seeing.beyond.target/tools/needle/db_human/protozoa/protozoa.ncbi.february.3.2018.fasta /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.fa
[main] Real time: 0.226 sec; CPU: 0.229 sec
Open look up table VIPR.table.txt.gz
Results are here /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.SV.csv
Done!
Open look up table  fungi.table.txt
Results are here /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.SV.csv
Done!
Open look up table protozoa.table.txt
Results are here /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.SV.csv
Done!
-----------------------------------------------------
Map assembled contigs onto the entire TREE of life to verify specificity of assembled contigs
Warning: [blastn] Query is Empty!
Warning: [blastn] The parameter -max_target_seqs is ignored for output formats, 0,1,2,3. Use -num_descriptions and -num_alignments to control output
Warning: [blastn] Query is Empty!
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.BLAST.house.format.csv')
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.SV.csv')
Contigs which are only mapped to custom microbial database and is not mapped to TREE OF LIFE FROM BLAST
set([])
Contigs which are  mapped to custom microbial database and TREE OF LIFE FROM (BLAST)
set([])
0  filtered contigs are saved to /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.SV.filtered.csv
0 contigs which better match bacteria are detected. We report only ones which are different by at least 2% /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.SV.filtered.csv.putative.bacteria.csv
done!
Warning: [blastn] The parameter -max_target_seqs is ignored for output formats, 0,1,2,3. Use -num_descriptions and -num_alignments to control output
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.BLAST.house.format.csv')
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.SV.csv')
Contigs which are only mapped to custom microbial database and is not mapped to TREE OF LIFE FROM BLAST
set([])
Contigs which are  mapped to custom microbial database and TREE OF LIFE FROM (BLAST)
set(['k61_13'])
===================================
===================================
['k61_13', 'KY429839.1', 'Homo_sapiens_clone_NA12878_chr6_160219097_160219098_genomic_sequence', '96.698', '212.0', '232.0', '88.3619655172']
['k61_13', 'NC_007257.2', 'Leishmania_major_strain_Friedlin_complete_genome__chromosome_16', '161', '79', '232', '-103.797468354']
['k61_13', 'NC_007257.2', 'Leishmania_major_strain_Friedlin_complete_genome__chromosome_16', '161', '79', '232', '-103.797468354']
KY429839 NC_007257.2 Homo Leishmania -103.797468354 88.3619655172
genus---> Homo Leishmania
0  filtered contigs are saved to /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.SV.filtered.csv
0 contigs which better match bacteria are detected. We report only ones which are different by at least 2% /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.SV.filtered.csv.putative.bacteria.csv
done!
Warning: [blastn] The parameter -max_target_seqs is ignored for output formats, 0,1,2,3. Use -num_descriptions and -num_alignments to control output
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.BLAST.house.format.csv')
('Open ', '/u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.SV.csv')
Contigs which are only mapped to custom microbial database and is not mapped to TREE OF LIFE FROM BLAST
set([])
Contigs which are  mapped to custom microbial database and TREE OF LIFE FROM (BLAST)
set(['k61_17', 'k61_47'])
===================================
===================================
['k61_17', 'NM_014491.4', 'Homo_sapiens_forkhead_box_P2_(FOXP2)_transcript_variant_1_mRNA', '100.0', '123.0', '203.0', '60.5911330049']
['k61_17', 'NC_006038.1', 'Kluyveromyces_lactis_strain_NRRL_Y-1140_chromosome_B_complete_sequence', '166', '38', '203', '-336.842105263']
['k61_17', 'NC_006038.1', 'Kluyveromyces_lactis_strain_NRRL_Y-1140_chromosome_B_complete_sequence', '166', '38', '203', '-336.842105263']
NM_014491 NC_006038.1 Homo Kluyveromyces -336.842105263 60.5911330049
genus---> Homo Kluyveromyces
===================================
===================================
['k61_47', 'KY429839.1', 'Homo_sapiens_clone_NA12878_chr6_160219097_160219098_genomic_sequence', '97.021', '235.0', '232.0', '98.2755818966']
['k61_47', 'NC_006037.1', 'Kluyveromyces_lactis_strain_NRRL_Y-1140_chromosome_A_complete_sequence', '195', '38', '232', '-413.157894737']
['k61_47', 'NC_006037.1', 'Kluyveromyces_lactis_strain_NRRL_Y-1140_chromosome_A_complete_sequence', '195', '38', '232', '-413.157894737']
KY429839 NC_006037.1 Homo Kluyveromyces -413.157894737 98.2755818966
genus---> Homo Kluyveromyces
0  filtered contigs are saved to /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.SV.filtered.csv
0 contigs which better match bacteria are detected. We report only ones which are different by at least 2% /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.SV.filtered.csv.putative.bacteria.csv
done!
Contigs after filtering are here
1 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.virus.megahit.contigs.SV.filtered.csv
1 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.fungi.megahit.contigs.SV.filtered.csv
1 /u/project/zarlab/serghei/op_test/test/needle.out/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.unmapped.fastq.protozoa.megahit.contigs.SV.filtered.csv
Success!!!
------------- (4) rDNA dosage -------------
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 10569 reads
prepare , /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.cat.rDNA.fastq
map reads to , /u/project/zarlab/serghei/op_test/test/LUAD_Mex_P80_T-G3_L_152666.dedup.cleaned.sort.rDNA.bam
(ERR): "/u/project/zarlab/serghei/code/seeing.beyond.target/db.human/rDNA.fasta" does not exist or is not a Bowtie 2 index
Exiting now ...
[warning] samtools mpileup option `u` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[mpileup] 1 samples in 1 input files
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
Traceback (most recent call last):
  File "/u/project/zarlab/serghei/code/seeing.beyond.target/bcf2maf.py", line 90, in <module>
    average_diversity=sum(sdi_list) / float(len(sdi_list))
ZeroDivisionError: float division by zero
------------- (5) MT dosage and diversity -------------
(ERR): "/u/project/zarlab/serghei/code/seeing.beyond.target/db.human/MT.fasta" does not exist or is not a Bowtie 2 index
Exiting now ...
[warning] samtools mpileup option `u` is functional, but deprecated. Please switch to using bcftools mpileup in future.
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
[mpileup] 1 samples in 1 input files
Traceback (most recent call last):
  File "/u/project/zarlab/serghei/code/seeing.beyond.target/bcf2maf.py", line 90, in <module>
    average_diversity=sum(sdi_list) / float(len(sdi_list))
ZeroDivisionError: float division by zero
------------- (optional) TE elements. Only for OncoPanel data -------------
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 7533320 reads
7533320 reads; of these:
  7533320 (100.00%) were unpaired; of these:
    7092099 (94.14%) aligned 0 times
    149499 (1.98%) aligned exactly 1 time
    291722 (3.87%) aligned >1 times
5.86% overall alignment rate
------------- (6) T and B cell repetoire profiling -------------
[W::sam_hdr_parse] Duplicated sequence 'M21626|TRAV14/DV4*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV23/DV6*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV29/DV5*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV36/DV7*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000661|TRAV38-2/DV8*01|Homo'
2313471 reads; of these:
  2313471 (100.00%) were unpaired; of these:
    2313431 (100.00%) aligned 0 times
    24 (0.00%) aligned exactly 1 time
    16 (0.00%) aligned >1 times
0.00% overall alignment rate
[W::sam_hdr_parse] Duplicated sequence 'M21626|TRAV14/DV4*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV23/DV6*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV29/DV5*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV36/DV7*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000661|TRAV38-2/DV8*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'M21626|TRAV14/DV4*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV23/DV6*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV29/DV5*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000660|TRAV36/DV7*01|Homo'
[W::sam_hdr_parse] Duplicated sequence 'AE000661|TRAV38-2/DV8*01|Homo'
[warning] samtools mpileup option `u` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[mpileup] 1 samples in 1 input files
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
------------- (7) Calculate off target coverage -------------






/u/home/s/serghei/project/code/seeing.beyond.target/sbt.sh: line 331: syntax error near unexpected token `else'
/u/home/s/serghei/project/code/seeing.beyond.target/sbt.sh: line 331: `else'
```
