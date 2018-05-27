import pysam
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('inbam', help='Mapped reads in bam format')
ap.add_argument('out_fastq', help='fastq with reads from chr1 228743523, 228781906')
args = ap.parse_args()


samfile = pysam.AlignmentFile(args.inbam, "rb")
out=open(args.out_fastq,"w")




for read in samfile.fetch('1', 228743523, 228781906):
	out.write("@"+read.query_name)
	out.write("\n")
	out.write(read.query_sequence)
	out.write("\n")
	out.write("+")
	out.write(read.qual)
	out.write("+")


samfile.close()

