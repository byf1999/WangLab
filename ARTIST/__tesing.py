import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
df = pd.read_excel(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM\insertion_distribution.xlsx')
plt.scatter(df['log'], list(range(len(df))))
df1 = df['Unnamed: 12']
df1 = df1[~np.isnan(df1)]

ax = plt.subplot()
ax.bar(list(range(len(df1))),df1)
ax.set_xlabel('Log2(# of insertions)',size=12)
ax.set_ylabel('# of genes',size=12)
ax.set_xticks(list(range(len(df1))))
plt.show()
ax.get_figure().savefig(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM\pics\distribution.svg')

with open(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\GSEA\gene_set1.txt', 'w') as fo:
    for tpe in set(df['eggNOG_class_annotation']):
        genes = df['#GeneID'][df['eggNOG_class_annotation'] == tpe]
        line = '\t'.join(genes)
        fo.write(f'{tpe}\tna\t{line}\n')