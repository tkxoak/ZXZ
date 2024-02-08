import pandas as pd
from collections import Counter

data = pd.read_csv('./GSE132465_annotation.txt', sep='\t', header=None)

word_counts = Counter(data[5])

with open('cell_counts.txt', 'w') as f:
    for word, count in word_counts.items():
        f.write(word + '\t' + str(count) + '\n')