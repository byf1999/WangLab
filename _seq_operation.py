from Bio import SeqIO
import pandas as pd
import os


def _extract_sequence(gene_list, annotation_file):
    '''extract sequences from a fasta file in the gene list'''

    def open_excel(file_path):
        return list(pd.read_excel(file_path, header=None, dtype=str).iloc[:, 0])

    def open_txt(file_path):
        with open(file_path, 'r') as input_:
            return [line.strip('\n') for line in input_.readlines()]

    operation_dict = {'xlsx': open_excel, 'xls': open_excel,
                      'txt': open_txt}

    suffix = os.path.splitext(gene_list)[-1][1:]
    gene_list = operation_dict.get(suffix, KeyError)(gene_list)

    sequence_list = [(record.id, record.seq) for record in
                     SeqIO.parse(annotation_file, os.path.splitext(annotation_file)[-1][1:])
                     if record.id in gene_list]

    return sequence_list


extract_sequence = _extract_sequence

if __name__ == '__main__':
    out = extract_sequence(r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\test_data\gene_list.xlsx',
                           r'F:\LearningFiles\Master\18.LearningMaterial\Python\wanglab\test_data\annotation.fasta')

    fo_path = r''
    with open(fo_path) as fo:
        for n, seq in out:
            fo.write('>' + n + '\n')
            fo.write(str(seq) + '\n')
