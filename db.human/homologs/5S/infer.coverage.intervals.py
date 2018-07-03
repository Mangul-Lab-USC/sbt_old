import csv
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('file', help='sorted bam file with mapped reads')
ap.add_argument('chr', help='chromosome number')
ap.add_argument('out', help='out')

args = ap.parse_args()





file=open(args.file)
reader=csv.reader(file)
pos=0



coordinates=[]



for line in reader:
	pos_previous=pos
	pos=int(line[1])
	if pos-pos_previous>5:
		coordinates.append((pos,pos_previous))


coordinates.append((pos,0))


fileOut=open(args.out,"w")

	
for i in range(1,len(coordinates)-1):
	fileOut.write(args.chr+","+str(coordinates[i-1][0])+","+str(coordinates[i][1]))
	fileOut.write("\n")

fileOut.close()



