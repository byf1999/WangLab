'''kegg API by Wang lab
URL form:
https://rest.kegg.jp/<operation>/<argument>[/<argument2[/<argument3> ...]]
<operation> = info | list | find | get | conv | link | ddi

info – display database release information and linked db information
list – obtain a list of entry identifiers and associated definition
find – find entries with matching query keyword or other query data
get – retrieve given database entries
conv – convert KEGG identifiers to/from outside identifiers
link – find related entries by using database cross-references
ddi – find adverse drug-drug interactions
'''
from Bio.KEGG import REST
import re
import time


# operation = ['info', 'list', 'find', 'get', 'conv', 'link', 'ddi']
# PREFIX = 'https://rest.kegg.jp/'


def _get_gene_name(gene_name, org=None):
    '''if more than one genes are found, return -1, which means the description should be more specific.
    if the arg org is defined, the organism will be checked.
    the gene_name will return only if there is only one item found in the specified organism.

    :return a tuple of gene name and organism
    '''
    gene_names = REST.kegg_find('genes', gene_name).read().rstrip().split("\n")
    lines = [line for line in gene_names if line.split('\t')[0].startswith(org)] if org else gene_names
    if len(lines) == 1:
        gene_name = lines[0].split('\t')[0]
        org = org if org else gene_name.split(':')[0]
        return gene_name
    time.sleep(0.1)
    return None  # Gene not found or more than one genes are found.


def _get_pathway_name(pathway_name):
    pathways = REST.kegg_find('pathway', pathway_name).read().rstrip().split('\n')
    time.sleep(0.1)
    return pathways[0].split('\t')[1]


def _gene_name_to_ko(gene_name_list, org=None):
    '''transfer a list of gene name to ko list'''
    gene_names = [_get_gene_name(gene_name, org=org) for gene_name in gene_name_list]
    kos = []
    for gene_name in gene_names:
        if gene_name:
            ko_links = REST.kegg_link('ko', gene_name).read().rstrip().split('\n')
            time.sleep(0.1)
            kos.append('; '.join([line.split('\t')[1] for line in ko_links if line]))
        else:
            kos.append('')
    return kos


def _gene_name_to_pathway(gene_name_list, org=None):
    '''transfer a list of gene name to pathway list'''
    gene_names = [_get_gene_name(gene_name, org=org) for gene_name in gene_name_list]
    pathways = []
    for gene_name in gene_names:
        if gene_name:
            pathway_link = REST.kegg_link('pathway', gene_name).read().rstrip().split('\n')
            time.sleep(0.1)
            pathways.append('; '.join([line.split('\t')[1] for line in pathway_link if line]))
        else:
            pathways.append('')
    return pathways

def _ko_to_pathway_with_org(ko_list: list, org: str) -> list:
    '''transfer a list of ko to pathway list'''
    l = len(ko_list)
    ko_org_link = REST.kegg_link(org, 'ko').read().rstrip().split('\n')
    org_path_link = REST.kegg_link('pathway', org).read().rstrip().split('\n')
    index_dict = {}
    for ind, ko in enumerate(ko_list):
        if ko in index_dict:
            index_dict[ko].append(ind)
        else:
            index_dict[ko] = [ind]

    gene_list = [None] * l
    for line in ko_org_link:
        ko, gene = line.split('\t')
        if ko[3:] in ko_list:
            ind = ko_list.index(ko[3:])
            gene_list[ind] = gene

    pathway_list = [[]] * l
    for line in org_path_link:
        gene, path = line.split('\t')
        if gene in gene_list:
            ind = gene_list.index(gene)
            pathway_list[ind] = pathway_list[ind] + [path[5:]]

    for ko, index in index_dict.items():
        if len(index) > 1:
            dta = pathway_list[ko_list.index(ko)]
            for ind in index[1:]:
                pathway_list[ind] = dta

    return pathway_list

def _ko_to_pathway_without_org(ko_list: list) -> list:
    '''transfer a list of ko to pathway list'''
    l = len(ko_list)
    index_dict = {}
    for ind, ko in enumerate(ko_list):
        if ko in index_dict:
            index_dict[ko].append(ind)
        else:
            index_dict[ko] = [ind]

    pathway_link = REST.kegg_link('pathway', 'ko').read().rstrip().split('\n')
    pathway_list = [[]] * l
    for line in pathway_link:
        ko, path = line.split('\t')
        if path.startswith('path:map') and ko[3:] in ko_list:
            ind = ko_list.index(ko[3:])
            pathway_list[ind] = pathway_list[ind] + [path[5:]]
    for ko, index in index_dict.items():
        if len(index) > 1:
            dta = pathway_list[ko_list.index(ko)]
            for ind in index[1:]:
                pathway_list[ind] = dta
    return pathway_list


def _pathway_to_description(pathway_list):
    descriptions = [_get_pathway_name(re.search('[a-zA-Z:]*(\d*)', pathway).group(1)) for pathway in pathway_list]
    return descriptions


gene_name_to_pathway = _gene_name_to_pathway
gene_name_to_ko = _gene_name_to_ko
pathway_to_description = _pathway_to_description
def ko_to_pathway(ko_list: list, org=None):
    if not org:
        return _ko_to_pathway_without_org(ko_list)
    elif org:
        return _ko_to_pathway_with_org(ko_list, org)
