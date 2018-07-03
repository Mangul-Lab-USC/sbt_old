. /u/local/Modules/default/init/modules.sh
module load jellyfish 
jellyfish count -m 75 -s 100M -t 10 -C HSU13369.fasta 
jellyfish dump mer_counts.jf | grep -v ">" | awk '{i+=1;print ">"i"\n"$1 }' >HSU13369.kmers.75bp.fasta
 
