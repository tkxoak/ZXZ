#
# with open('./CIBERSORT/GSE132465_6_bulk.txt', 'r') as input_file:
#
#     lines = [line.strip().split('\t') for line in input_file.readlines()]
#
#
# transposed_lines = list(map(list, zip(*lines)))



with open('/home/xianzhen/R/monocle/GSE132465_count_matrix.txt', 'r') as file:
    lines = file.readlines()

lines = [line.strip().split('\t') for line in lines]

max_columns = max(len(line) for line in lines)

for i in range(len(lines)):
    if len(lines[i]) < max_columns:
        lines[i].extend([''] * (max_columns - len(lines[i])))

transposed = list(map('\t'.join, zip(*lines)))

with open('/home/xianzhen/R/monocle/GSE132465_count.txt', 'w') as file:
    file.write('\n'.join(transposed))
