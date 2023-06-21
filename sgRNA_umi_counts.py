# usage: python sgRNA_umi_counts.py read1 read2 out_prefix
from sys import argv
import gzip
from Bio.SeqIO import parse
import pandas as pd

read1=argv[1]#'tmp_R1.fastq.gz'
read2=argv[2]#'tmp_R2.fastq.gz'
out=argv[3] # out prefix

# read fasta, grna seq:id mapping
grna={}
for rec in parse('data/sgRNA_20nt.fa','fasta'):
    grna[str(rec.seq)]=str(rec.id)
   
n=0 # line number
r=0 # read number

## process read 2, store as read no:grna id
grnas={}
for line in gzip.open(read2,'rt'):
    n+=1
    if n%4==2:
        r+=1
        for s in grna:
            if s in line:
                grnas[r]=[grna[s]]
                break
                
## process read1                                
# 3 cell barcodes
cb=pd.read_csv('data/cb.csv')

# select barcode based on hd
def barcode(bc,b): # string, list of barcodes
    if bc in b:
        return bc
    else:
        for s in b:
            n=0
            for i in range(len(s)): # hamming dist
                if s[i]!=bc[i]:
                    n+=1
            if n==1: # hd=1
                return s
                break

n=0 # line number
r=0 # read number

# barcodes dict, read no: [sgrna,barcode,umi]
for line in gzip.open(read1,'rt'):
    n+=1
    if n%4==2:
        r+=1
        if r in grnas: # reads with grna
            # read barcodes
            bc1= barcode(line[0:9],   cb.b1.tolist()) 
            bc2= barcode(line[21:30], cb.b2.tolist()) 
            bc3= barcode(line[43:52], cb.b3.tolist()) 
            
            # combine barcodes
            if bc1 !=None and bc2 !=None and bc3 !=None: # skip missing barcodes
                bc='_'.join([bc1,bc2,bc3])
                grnas[r]+=[bc,line[52:60]]

                
## merge reads
# remove reads with missing cell barcodes
a=pd.DataFrame.from_dict(grnas,orient='index', columns=['grna','barcode','umi'])
n1=len(a)

# remove umi duplicates: count unique cell barcode, gRNA, umi
a=a.value_counts().reset_index()
n2=len(a)

print(read2,'UMI ratio: %.2f'%(n1/n2))

# sum umi counts
a=a[['barcode','grna']].value_counts().reset_index()

# save data
a.columns=['barcode','sgRNA','counts']
a.to_csv(out+'_umi_counts.csv.gz',index=False)

