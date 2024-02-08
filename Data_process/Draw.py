

# with open('./count.txt', 'r') as f_in:
#     header = f_in.readline().strip().split('\t')  # read file header
#     smc_cols = [i for i, col in enumerate(header) if col.startswith("SMC06-T")  or col.startswith("SMC07-T")]
#     smccol = [0]
#     smccol.extend(smc_cols)
#     with open('./CIBERSORT/GSE132465_ref.txt', 'w') as f_out:
#         f_out.write('\t'.join(header[i] for i in smccol) + '\n')
#         for line in f_in:
#             cols = line.strip().split('\t')
#             smc_data = [cols[i] for i in smccol]  # draw SMC01-T data
#             f_out.write('\t'.join(smc_data) + '\n')  #


with open('./count.txt', 'r') as f_in:
    header = f_in.readline().strip().split('\t')  # read file header
    smc_cols = [i for i, col in enumerate(header) if col.startswith("SMC01-T")  or col.startswith("SMC03-T") or col.startswith("SMC04-T") or col.startswith("SMC06-T") or col.startswith("SMC07-T")
    or col.startswith("SMC10-T")]
    smccol = [0]
    smccol.extend(smc_cols)
    with open('/home/xianzhen/PycharmProjects/dnad/MSI/GSE132465_6_ref.txt', 'w') as f_out:
        f_out.write('\t'.join(header[i] for i in smccol) + '\n')
        for line in f_in:
            cols = line.strip().split('\t')
            smc_data = [cols[i] for i in smccol]  # draw SMC01-T data
            f_out.write('\t'.join(smc_data) + '\n')  #