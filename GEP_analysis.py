import pandas as pd
import numpy as np
import torch

from SAE.utils import counts2TPM, ProcessInputData
from SAE.train import train_model, predict
import seaborn as sns
import matplotlib.pyplot as plt
from SAE.simulation import generate_simulated_data
from sklearn.preprocessing import MinMaxScaler
from scipy import stats
import matplotlib.colors as colors
import anndata

def CCCscore(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    covariance = np.cov(x, y)[0][1]
    x_var = np.var(x)
    y_var = np.var(y)
    ccc = (2 * covariance) / (x_var + y_var + (x_mean - y_mean) ** 2)
    return ccc

test_data = pd.read_csv('~/zxz/Dnad/GSE178341_limma/GSE178341_MSI-H_test_3.txt',sep='\t',index_col=0)
data8k = pd.read_csv('~/zxz/Dnad/MSI/GSE132465_6_ref_2.txt',sep='\t',index_col=0)#
#data8k = counts2TPM(data8k,genelen='~/PycharmProjects/SAE-main/data/GeneLength.txt')

original_sigm = data8k.groupby('Cell_type').mean()



# crr_dir = os.getcwd()       # get the cwd
# os.chdir(crr_dir + '/src')

original_sigm = np.log2(original_sigm+1)
train_x, train_y, test_x, genename, celltypes, samplename = \
    ProcessInputData('/home/zxz/zxz/SAE-main/GSE132465_6_ref.h5ad', test_data, sep='\t',datatype='counts')
#model = train_model(train_x,train_y,batch_size=128,model_name= "GEP_als")
model = torch.load('GEP_als.pth')
train_sigm = model.sigmatrix().cpu().detach().numpy()
Sigm, Pred = \
    predict(test_x=test_x,genename=genename,celltypes=celltypes,samplename=samplename,
            model=model,
            adaptive=True,mode='overall')
train_sigm = pd.DataFrame(model.sigmatrix().cpu().detach().numpy(),columns=genename, index=celltypes)
inter = Sigm.columns.intersection(original_sigm.columns)
Sigm = Sigm[inter]
original_sigm = original_sigm[inter]
train_sigm = train_sigm[inter]
inter = Sigm.index.intersection(original_sigm.index)
Sigm = Sigm.loc[inter]
original_sigm = original_sigm.loc[inter]
train_sigm= train_sigm.loc[inter]
mms = MinMaxScaler()
Sigm = mms.fit_transform(Sigm.T).T
original_sigm = mms.fit_transform(original_sigm.T).T
train_sigm = mms.fit_transform(train_sigm.T).T

# sns.scatterplot(original_sigm[5,:], Sigm[5,:],marker='.')
# plt.show()
# sns.scatterplot(x=original_sigm[4,:], y=Sigm[4,:])#,marker='.'
# ccc4 = CCCscore(original_sigm[4,:], Sigm[4,:])
# print("ccc4=" + str(ccc4))
# plt.show()
# sns.scatterplot(x=original_sigm[3,:], y=Sigm[3,:])
# ccc3 = CCCscore(original_sigm[3,:], Sigm[3,:])
# print("ccc3=" + str(ccc3))
# plt.show()
# sns.scatterplot(x=original_sigm[2,:], y=Sigm[2,:])
# ccc2 = CCCscore(original_sigm[2,:], Sigm[2,:])
# print("ccc2=" + str(ccc2))
# plt.show()
# sns.scatterplot(x=original_sigm[1,:], y=Sigm[1,:])
# ccc1 = CCCscore(original_sigm[1,:], Sigm[1,:])
# print("ccc1=" + str(ccc1))
# plt.show()
# sns.scatterplot(x=original_sigm[0,:], y=Sigm[0,:])
# ccc0 = CCCscore(original_sigm[0,:], Sigm[0,:])
# print("ccc0=" + str(ccc0))
# plt.show()

# sns.scatterplot(x=original_sigm[4,:], y=train_sigm[4,:],marker='.')
# plt.show()
# sns.scatterplot(x=original_sigm[3,:], y=train_sigm[3,:],marker='.')
# plt.show()
# sns.scatterplot(x=original_sigm[2,:], y=train_sigm[2,:],marker='.')
# plt.show()
# sns.scatterplot(x=original_sigm[1,:], y=train_sigm[1,:],marker='.')
# plt.show()
# sns.scatterplot(x=original_sigm[0,:], y=train_sigm[0,:],marker='.')
# plt.show()

plt.scatter(x=original_sigm[5,:], y=Sigm[5,:],s=3)#,marker='.'
plt.title('B cells')
ccc4 = CCCscore(original_sigm[5,:], Sigm[5,:])
print("ccc5=" + str(ccc4))
plt.show()
plt.scatter(x=original_sigm[4,:], y=Sigm[4,:],s=3)
plt.title('Epithelial cells')
ccc3 = CCCscore(original_sigm[4,:], Sigm[4,:])
print("ccc4=" + str(ccc3))
plt.show()
plt.scatter(x=original_sigm[3,:], y=Sigm[3,:],s=3)
plt.title('Mast cells')
ccc3 = CCCscore(original_sigm[3,:], Sigm[3,:])
print("ccc3=" + str(ccc3))
plt.show()
plt.scatter(x=original_sigm[2,:], y=Sigm[2,:],s=3)
plt.title('Myeloids')
ccc2 = CCCscore(original_sigm[2,:], Sigm[2,:])
print("ccc2=" + str(ccc2))
plt.show()
plt.scatter(x=original_sigm[1,:], y=Sigm[1,:],s=3)
plt.title('Stromal cells')
ccc1 = CCCscore(original_sigm[1,:], Sigm[1,:])
print("ccc1=" + str(ccc1))
plt.show()
plt.scatter(x=original_sigm[0,:], y=Sigm[0,:],s=3)
plt.title('T cells')
ccc0 = CCCscore(original_sigm[0,:], Sigm[0,:])
print("ccc0=" + str(ccc0))
plt.show()