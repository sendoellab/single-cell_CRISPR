import scanpy as sc
import pandas as pd

a=sc.read_h5ad('../seurat_new/out/wta11to18_filtered_seurat_merge_raw.h5ad')

obs=pd.DataFrame()

# estimate doublet for each sample separately
obs=pd.DataFrame()
for sample in 'wta11 wta12 wta13 wta14 wta15 wta16 wta17 wta18'.split(): #
    print(sample)
    wta = a[a.obs['sample'] == sample]
    sc.pp.filter_genes(wta, min_cells=5)
    sc.external.pp.scrublet(wta, expected_doublet_rate=0.2, n_prin_comps=40)
    sc.external.pl.scrublet_score_distribution(wta, save='_'+sample+'_pc40.pdf')
    obs=pd.concat([wta.obs,obs])

# add sgRNA annotations for cells
pa='gRNA_counts' 

an11=pd.read_csv(pa+'WTA11_S1_umi_counts_anno.csv.gz')
an11['cell']=an11.id.astype(str)+'_1'

an12=pd.read_csv(pa+'WTA12_S2_umi_counts_anno.csv.gz')
an12['cell']=an12.id.astype(str)+'_2'

an13=pd.read_csv(pa+'WTA13_S3_umi_counts_anno.csv.gz')
an13['cell']=an13.id.astype(str)+'_3'

an14=pd.read_csv(pa+'WTA14_S4_umi_counts_anno.csv.gz')
an14['cell']=an14.id.astype(str)+'_4'

an15=pd.read_csv(pa+'WTA15_S5_umi_counts_anno.csv.gz')
an15['cell']=an15.id.astype(str)+'_5'

an16=pd.read_csv(pa+'WTA16_S7_umi_counts_anno.csv.gz')
an16['cell']=an16.id.astype(str)+'_6'

an17=pd.read_csv(pa+'WTA17_S8_umi_counts_anno.csv.gz')
an17['cell']=an17.id.astype(str)+'_7'

an18=pd.read_csv(pa1+'/WTA18_S6_umi_counts_anno.csv.gz')
an18['cell']=an18.id.astype(str)+'_8'

an=pd.concat([an11, an12, an13, an14, an15, an16, an17, an18]).set_index('cell')
an['gRNA']=an.sgRNA.apply(lambda x:x.split('_')[0]) # add grna anno

# combine doublet info with sgRNA anno
obs1=obs.join(an)

# select sgRNA based on quantile 0.99 cutoff and doublet filter
obs2=obs1[(obs1.counts>obs1.q99) | (obs1.sgrna_detected==1)]
obs2=obs2[obs2.predicted_doublet==False]

obs2.iloc[:,4:].drop('barcode',axis=1).to_csv('out/wta11to18_filtered_seurat_merge_raw_scrublet_pc40_grna_anno_q99_nodoublet.csv')