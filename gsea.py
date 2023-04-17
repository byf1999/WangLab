import pandas as pd
import os
from _kegg import *
from Bio.KEGG import REST

def make_kegg_gene_set(fi_path, fo_path=None):
    '''given an expression file and output a .gmt file used in GSEA'''
    df = pd.read_table(fi_path, sep='\t')
    ko_list = list(set(df.iloc[:,1]))
    transform_dict = {}
    for name, ko in zip(df.iloc[:,0], df.iloc[:,1]):
        if ko in transform_dict:
            transform_dict[ko] += f'\t{name}'
        else:
            transform_dict[ko] = name


    pathway_list = ko_to_pathway(ko_list)
    all_pathways = set([item for l in pathway_list for item in l])

    if not fo_path:
        fo_path = os.path.join(os.path.dirname(fi_path), 'kegg_gene_set.gmt')
    with open(fo_path, 'w') as fo:
        for t in all_pathways:
            pathway_link = REST.kegg_link('ko', t).read().rstrip().split('\n')
            list_of_pathway = [line.split('\t')[1][3:] for line in pathway_link]
            intersect = set(ko_list).intersection(list_of_pathway)

            line = '\t'.join([transform_dict[map] for map in intersect])
            fo.write(f'{t}\tna\t{line}\n')
if __name__ == '__main__':
    fi_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\GSEA\expression1.txt'
    make_kegg_gene_set(fi_path)