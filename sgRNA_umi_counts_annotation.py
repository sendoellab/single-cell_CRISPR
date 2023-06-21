# usage: python 
from glob import glob
import pandas as pd

#f='WTA10_S6'
for f in glob('*_umi_counts.csv.gz'):

    a=pd.read_csv(f)
    # cells with number of gRNA
    sgrna_detected=a.barcode.value_counts().to_frame('sgrna_detected')

    # select gRNA [['gRNA','barcode','counts']]
    aa=a.groupby(['barcode','sgRNA']).sum().sort_values(by='counts',ascending=False)

    # calculate quantile
    aq=aa.groupby(['barcode']).quantile(0.99)
    aq.columns=['q99']

    # select top gRNA
    ab=aa.groupby(level=0).head(1).reset_index().set_index('barcode')
    ab=ab.join(aq)

    # expression filter
    #ab=ab[ab.counts>ab.q95]

    # cell barcodes mapping
    m=pd.read_table('data/barcode_mapping.tsv.gz',index_col=0)#.astype(int)
    ab=ab.join(m).dropna()
    ab['id']=ab['id'].astype(int)

    # save annotation
    f1=f.split('/')[-1].replace('_counts','_counts_anno')
    ab.join(sgrna_detected).round(3).to_csv(f1)