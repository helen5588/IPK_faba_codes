
#python heatmap.raw.py 1_607075117-809001340.1k.geno.vcf.012 All.mat.tree.id 
#1_607075117-809001340.1k.geno.vcf.012 can be produced by vcftools --012
#All.mat.tree.id:samplename\tclassification
#Arrisot_1_EP1-1-1-1     spring
import sys
import seaborn as sns
import pandas as pd 
import numpy as np 
from pandas import Series,DataFrame
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
train = pd.read_table(sys.argv[1],delimiter="\t",index_col=0,header=None)
bb=pd.read_table(sys.argv[1]+".indv",delimiter="\t",index_col=0,header=None)
train.index = bb.index
cc=pd.read_table(sys.argv[2],delimiter="\t",index_col=0,header=None)
dd=pd.read_table(sys.argv[1]+".pos",delimiter="\t",header=None)
train.columns = dd[1]
tra = train.reindex(cc.index)
colors = ["#992d22","#c27c0e","#27ae60","#1f618d"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12), gridspec_kw={'width_ratios': [0.02, 1],'wspace': 0.01})
unique_values = sorted(set(cc[1]))
colors2 = sns.color_palette("coolwarm", n_colors=len(cc[1].unique()))#RdBu
row_colors = cc[1].map(dict(zip(cc[1].unique(), colors2)))
a = sns.clustermap(tra,cmap=sns.color_palette(colors),annot=False, row_cluster=True,col_cluster=False, yticklabels=False,row_colors=row_colors, cbar_kws={},cbar_pos=None)

row_order = a.dendrogram_row.reordered_ind
row_names = tra.index[row_order]

a.ax_heatmap.set_xlabel('')
a.ax_heatmap.set_ylabel('')
plt.subplots_adjust(wspace=0)
plt.savefig(sys.argv[1]+".png")
plt.show()
