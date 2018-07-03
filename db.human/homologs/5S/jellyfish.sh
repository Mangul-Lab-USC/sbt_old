. /u/local/Modules/default/init/modules.sh
module load jellyfish 
jellyfish count -m 75 -s 100M -t 10 -C 5S.fasta 
jellyfish dump mer_counts.jf | grep -v ">" | awk '{i+=1;print ">"i"\n"$1 }' >5S.kmers.75bp.fasta
 
