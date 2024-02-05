import os
from SAE.simulation import generate_simulated_data
from SAE.utils import ProcessInputData
import csv
import random
import torch
from SAE.train import adaptive_stage
import numpy as np
from tqdm import tqdm
import pandas as pd
from SAE.model import  device
def shuffle_tab_separated_file(file_path):
    with open(file_path, 'r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        lines = list(reader)
    header = lines[0]
    data_lines = lines[1:]
    random.shuffle(data_lines)
    shuffled_lines = [header] + data_lines
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(shuffled_lines)


def test_predict(test_x, genename, celltypes, samplename,
            model_name=None, model=None,
            adaptive=True, mode='overall'):
    if adaptive is True:
        if mode == 'high-resolution':
            TestSigmList = np.zeros((test_x.shape[0], len(celltypes), len(genename)))
            TestPred = np.zeros((test_x.shape[0], len(celltypes)))
            print('Start adaptive training at high-resolution')
            for i in tqdm(range(len(test_x))):
                x = test_x[i, :].reshape(1, -1)
                if model_name is not None and model is None:
                    model = torch.load(model_name + ".pth")
                elif model is not None and model_name is None:
                    model = torch.load("model.pth")
                decoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'decoder' in n]}]
                encoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'encoder' in n]}]
                optimizerD = torch.optim.Adam(decoder_parameters, lr=1e-4)
                optimizerE = torch.optim.Adam(encoder_parameters, lr=1e-4)
                test_sigm, loss, test_pred = adaptive_stage(model, x, optimizerD, optimizerE, step=300, max_iter=3)
                # mms = MinMaxScaler()#
                TestSigmList[i, :, :] = test_sigm
                TestPred[i, :] = test_pred
            TestPred = pd.DataFrame(TestPred, columns=celltypes, index=samplename)
            CellTypeSigm = {}
            for i in range(len(celltypes)):
                cellname = celltypes[i]
                sigm = TestSigmList[:, i, :]
                # with open('scaler.pkl', 'rb') as f:
                #     scaler = pickle.load(f)
                # sigm = scaler.inverse_transform(sigm.T)  #
                # sigm = sigm.T  #
                sigm = pd.DataFrame(sigm, columns=genename, index=samplename)
                CellTypeSigm[cellname] = sigm
            print('Adaptive stage is done')

            return CellTypeSigm, TestPred

        elif mode == 'overall':
            if model_name is not None and model is None:
                try:
                    model = torch.load(model_name + ".pth")
                    print("Model loaded successfully.")
                except Exception as e:
                    print(f"Error loading the model: {e}")
            elif model is not None and model_name is None:
                try:
                    model = torch.load("/home/zxz/zxz/SAE-main/Test/model.pth")
                    print("Model loaded successfully.")
                except Exception as e:
                    print(f"Error loading the model: {e}")
            decoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'decoder' in n]}]
            encoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'encoder' in n]}]
            optimizerD = torch.optim.Adam(decoder_parameters, lr=1e-4)
            optimizerE = torch.optim.Adam(encoder_parameters, lr=1e-4)
            print('Start adaptive training for all the samples')
            test_sigm, loss, test_pred = adaptive_stage(model, test_x, optimizerD, optimizerE, step=300, max_iter=3)
            print('Adaptive stage is done')
            test_sigm = pd.DataFrame(test_sigm, columns=genename, index=celltypes)
            test_pred = pd.DataFrame(test_pred, columns=celltypes, index=samplename)

            return test_sigm, test_pred

    else:
        if model_name is not None and model is None:
            model = torch.load(model_name + ".pth")
        elif model is not None and model_name is None:
            model = model
        print('Predict cell fractions without adaptive training')
        model.eval()
        model.state = 'test'
        data = torch.from_numpy(test_x).float().to(device)
        _, pred, _ = model(data)
        pred = pred.cpu().detach().numpy()
        pred = pd.DataFrame(pred, columns=celltypes, index=samplename)
        print('Prediction is done')
        return pred


def test_Deconvolution(simudata, real_bulk, sep='\t', variance_threshold=0.8,
                  datatype='counts', genelenfile=None,model = 'model.pth',
                  mode='overall', adaptive=True
                  ):

    train_x, train_y, test_x, genename, celltypes, samplename = \
        ProcessInputData(simudata, real_bulk, sep=sep, datatype=datatype, genelenfile=genelenfile,
                         variance_threshold=variance_threshold)
    print('Notice that you are using parameters: mode=' + str(mode) + ' and adaptive=' + str(adaptive))

    if adaptive is True:
        if mode == 'high-resolution':
            CellTypeSigm, Pred = \
                test_predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                        model_name = 'model',
                        adaptive=adaptive, mode=mode)
            return CellTypeSigm, Pred

        elif mode == 'overall':
            Sigm, Pred = \
                test_predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                        model_name = 'model',
                        adaptive=adaptive, mode=mode)
            return Sigm, Pred
    else:
        Pred = test_predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                       model_name = 'model',
                       adaptive=adaptive, mode=mode)
        Sigm = None
        return Sigm, Pred

folder_name = 'GSE119409'
os.makedirs(folder_name, exist_ok=True)

simdata = generate_simulated_data('/home/zxz/zxz/Dnad/MSI/GSE132465_6_ref_2.txt')
for i in range(1, 101):
    shuffle_tab_separated_file("/home/zxz/R/BULK_process/GSE119409_1.txt")
    file_name = f'Pred_{i}.txt'
    file_path = os.path.join(folder_name, file_name)
    Sigm, Pred = test_Deconvolution(simdata, "/home/zxz/R/BULK_process/GSE119409_1.txt", sep='\t',
                               datatype='counts', genelenfile='./GeneLength.txt',
                               mode='overall', adaptive=True
                               )
    Pred.to_csv(file_path, sep='\t', index=True)


