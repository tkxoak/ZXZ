from SAE import Deconvolution

Sigm, Pred = Deconvolution('~/R/GSE169147_2.txt', '~/R/GSE8056/GSE8056_1.txt',sep='\t',
                           datatype='counts', genelenfile='./GeneLength.txt',
                           mode='high-resolution', adaptive=True,
                           save_model_name = None)


for key, value in Sigm.items():
    value.to_csv(f'{key}.txt', index=False,sep='\t')

Pred.to_csv("Pred.txt", sep="\t", index=False)