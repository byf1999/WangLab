'''
https://rest.kegg.jp/link/etr/etr00010

'''

from _kegg import *
from _seq_operation import *
import pandas as pd
import numpy as np
import collections
import time

t0 = time.time()
df = pd.read_excel(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\gene_list.xlsx', header=0)
# find ko and pathway for each gene in the gene list
df['ko'] = gene_name_to_ko(df.iloc[:, 0], org='etr')
df['pathway'] = gene_name_to_pathway(df.iloc[:, 0], org='etr')

df.to_excel(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\output.xlsx', index=False)

# count genes for each pathway
all_pathways = []
for pathway in df['pathway']:
    if pathway:
        all_pathways += pathway.split('; ')

pathway_counter = collections.Counter(all_pathways)
pathway_list = [pathway for pathway, _ in pathway_counter.most_common()]
count_list = [count for _, count in pathway_counter.most_common()]
description_list = pathway_to_description(pathway_list)

df1 = pd.DataFrame(np.array([pathway_list, description_list, count_list]).transpose(),
                   columns=['pathway', 'description', 'count'])
df1.to_excel(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\pathways.xlsx', index=False)

print(fr'cost {time.time() - t0} s')

# extract_sequence
out = extract_sequence(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\test_data\gene_list.xlsx',
                       r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\test_data\annotation.fasta')

print(out)

df = pd.read_excel(r'F:\LearningFiles\Master\16.OtherPeople\DHY\user_ko (1).xlsx', header=None)
ko_list = list(df.iloc[:406, 1])
pathway_list = ko_to_pathway(ko_list, org='ppj')
pathway_list = ['; '.join(pathway) for pathway in pathway_list]
df2 = pd.DataFrame(np.array([ko_list, pathway_list]).transpose(),
                   columns=['ko', 'pathway'])
df2.to_excel(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\ko_pathways.xlsx', index=False)

out = []
for ko, pathways in zip(ko_list, pathway_list):
    for pathway in pathways:
        out.append([ko, pathway])

pd.DataFrame(np.array(out)).to_excel(r'')

# sam2wig
from tools import sam2wig

sam_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A.ve ARTIST scripts\A-V-Tn55_combined_R1_ART.SAM'
sam2wig(sam_file, chr_identifier='Contig00001')

# generate annotations
import pandas as pd
import numpy as np
import os
path = r'F:\LearningFiles\Master\19.A_ve_paper\Web_Report\Web_Report\customer_backup\function_annotation\01.general_database_annotation\01.combination_annotation'

df = pd.read_excel(os.path.join(path, 'Integrated_annotation.xlsx'), header=0)
n_col = len(df.columns)
new_df = []
l = 0
for i in range(int(df.iloc[-1, 0][-4:])):
    if df.iloc[l, 0].endswith(str(i+1)):
        new_df.append(list(df.iloc[l,]))
        l += 1
    else:
        new_df.append([f'GE{i+1:0>6}'] + ['' for _ in range(n_col - 1)])

new_df = pd.DataFrame(np.array(new_df), columns=df.columns)
new_df.to_excel(os.path.join(path, 'annotations.xlsx'), index=False)

# generate nlp dataset
'''
['Lysobacter_capsici_no_coli_homologs.fasta',
 'Lysobacter_capsici_test_data.fasta',
 'Lysobacter_capsici_train_data.fasta',
 'negative_test_data.fasta',
 'negative_train_data.fasta',
 'negative_Xanthomonas_data.fasta',
 'positive_test_data.fasta',
 'positive_train_data.fasta',
 'positive_Xantomonas_data.fasta']
'''
fi_path = r'C:\Users\byf1999\Downloads\T3ES_secretion_signal_analysis-main\T3ES_secretion_signal_analysis-main\data\100 N-terminal sequences of all datasets'
fo_path = r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\dl\dataset\nlp'
fi = 'negative_Xanthomonas_data.fasta'
count = 0
for record in SeqIO.parse(os.path.join(fi_path, fi), 'fasta'):
    fname = f'{count}.txt'
    count += 1
    seq = str(record.seq)
    with open(os.path.join(fo_path, 'test', 'xanthomonas',fi[:3], fname), 'w') as fo:
        fo.write(seq)