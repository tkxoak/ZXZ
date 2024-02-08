import anndata
import scanpy as sc
#bulk
adata = sc.read_h5ad("./GSE178341_Deseq2/GSE178341_MSS_test.h5ad")
gene_names = adata.var_names.tolist()
with open("./GSE178341_Deseq2/GSE178341_MSS_test_3.txt", "w") as f:
    f.write("\t" + "\t".join(gene_names) + "\n")
    for i, row in enumerate(adata.X):
        f.write(str(i + 1) + "\t" + "\t".join([str(x) for x in row.tolist()]) + "\n")

#celltype
# andata = anndata.read_h5ad("./CIBERSORT/CIBERSORT_GSE132465.h5ad")
#
# adata.obs.to_csv("./CIBERSORT/GSE132365_6_celltype.txt", sep="\t", index=False)