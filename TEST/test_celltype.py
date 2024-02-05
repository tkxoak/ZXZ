import os
import random
import anndata
import pandas as pd

from SAE import Deconvolution
import csv
import random

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



folder_name = 'GSE45404'
os.makedirs(folder_name, exist_ok=True)


for i in range(1, 101):
    shuffle_tab_separated_file("/home/zxz/R/BULK_process/GSE45404_1.txt")
    file_name = f'Pred_{i}.txt'
    file_path = os.path.join(folder_name, file_name)
    Sigm, Pred = Deconvolution('/home/zxz/zxz/Dnad/MSI/GSE132465_6_ref_2.txt', "/home/zxz/R/BULK_process/GSE45404_1.txt", sep='\t',
                               datatype='counts', genelenfile='./GeneLength.txt',
                               mode='overall', adaptive=True,
                               save_model_name = None)
    Pred.to_csv(file_path, sep='\t', index=True)



