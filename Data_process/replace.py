# # open file
# with open('./GSE178341_Deseq2/GSE178341_MSS_test_1.txt', 'r') as f1, open('/home/xianzhen/PycharmProjects/GSE178341/cluster.csv', 'r') as f2, open('./GSE178341_Deseq2/GSE178341_MSS_test_2.txt', 'w') as f3:
#     # Read the file and save the dictionary
#     replace_dict = {}
#     for line in f2:
#         elements = line.strip().split('\t')
#         replace_dict[elements[0]] = elements[2]
#
#     # replace
#     for line in f1:
#         elements = line.strip().split('\t')
#         if elements[0] in replace_dict:
#             elements[0] = replace_dict[elements[0]]
#         f3.write('\t'.join(elements) + '\n')

with open('/home/xianzhen/PycharmProjects/dnad/MSI/GSE132465_6_ref_1.txt', 'r') as f1, open('/home/xianzhen/PycharmProjects/dnad/GSE132465_annotation.txt', 'r') as f2, open('/home/xianzhen/PycharmProjects/dnad/MSI/GSE132465_6_ref_2.txt', 'w') as f3:
    # Read the file and save the dictionary
    replace_dict = {}
    for line in f2:
        elements = line.strip().split('\t')
        replace_dict[elements[0]] = elements[4]

    # replace
    for line in f1:
        elements = line.strip().split('\t')
        if elements[0] in replace_dict:
            elements[0] = replace_dict[elements[0]]
        f3.write('\t'.join(elements) + '\n')