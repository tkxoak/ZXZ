
import pandas as pd
import scanpy as sc

df1 = pd.read_csv('./gse132465_2TN_subtype.txt', sep='\t', index_col=0)
df2 = pd.read_csv('./gse132465_2TN_1_subtype.txt', sep='\t', index_col=0)

adata = sc.AnnData(df1)
adata = adata.concatenate(sc.AnnData(df2), join='outer', batch_key='batch')

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=5)

adata1 = adata[adata.obs['batch'] == '0']
adata2 = adata[adata.obs['batch'] == '1']

df_filtered1 = pd.DataFrame(adata1.X.toarray(), index=adata1.obs.index, columns=adata1.var.index)
df_filtered2 = pd.DataFrame(adata2.X.toarray(), index=adata2.obs.index, columns=adata2.var.index)

df_filtered1.to_csv('filtered_file1.txt', sep='\t')
df_filtered2.to_csv('filtered_file2.txt', sep='\t')

