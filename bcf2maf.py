from pysam import VariantFile
import argparse
import math
import scipy as sc
import os
from os.path import basename



def p(n, N):
    """ Relative abundance """
    if n is  0:
        return 0
    else:
        return (float(n)/N) * math.log(float(n)/N,2)

def sdi(data):
    if len(data)==0:
        return  "N/A"
    else:
        N = sum(data)
        return abs(-sum(p(n, N) for n in data if n is not 0))

ap = argparse.ArgumentParser()
ap.add_argument('input_bcf', help='--')
ap.add_argument('out', help='---')
args = ap.parse_args()

bcf_in = VariantFile(args.input_bcf)  # auto-detect input format


#test=[0.75,0.25]
#print "---->",sdi(test)  0.811278124459 which is CORRECT!



dict={}
pos=set()

for rec in bcf_in.fetch():
    #(0, 0, 8, 0)
    #VAF = (forward non-ref + reverse non-ref alleles) /  (forward ref alleles + reverse ref + forward non-ref + reverse non-ref alleles)
    
    

    ref=rec.info["DP4"][0]+rec.info["DP4"][1]
    non_ref=rec.info["DP4"][2]+rec.info["DP4"][3]
    if ref==0:
        VAF=1
    elif non_ref==0:
        VAF=0
    else:
        VAF=non_ref/ref

    print VAF

    if VAF!=0:
        if rec.pos not in pos:
            dict[rec.pos]=[]
            dict[rec.pos].append(VAF)
            pos.add(rec.pos)
        else:
            dict[rec.pos].append(VAF)
        






base=os.path.basename(args.input_bcf)
sample=os.path.splitext(base)[0]


sdi_list=[]

fileOut=open(args.out,"w")
fileOut.write("pos,diversity\n")
for key, value in dict.items():
    print value, sdi(value)

    fileOut.write(str(key)+","+str(sdi(value)))
    fileOut.write("\n")
    sdi_list.append(sdi(value))


fileOut.close()


average_diversity=sum(sdi_list) / float(len(sdi_list))


fileOut=open(args.out+".per.sample.csv","w")
fileOut.write("sample,average_diversity\n")
fileOut.write(sample+","+str(average_diversity) )
fileOut.write("\n")


print "done!"



