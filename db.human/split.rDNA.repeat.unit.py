from Bio import SeqIO


input_file='rDNA.fasta'

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
	if name=='gi|555853|gb|U13369.1|HSU13369':
		str_rDNA=sequence


#18S

input_file='18S.fasta'

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        str_18S=sequence


print str_rDNA
print "---"
print str_18S
print len(str_rDNA.split(str_18S)[0]), len(str_rDNA)

