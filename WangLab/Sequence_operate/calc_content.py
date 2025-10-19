import sys
from Bio import SeqIO


def read_file(file, fmt):
    seq_list = []
    if fmt == 't':
        with open(file, 'r') as fi:
            while True:
                content = fi.readline()
                if not content:
                    break
                if c := content.rstrip('\n'):
                    seq_list.append(c.split('\t'))
    elif fmt == ',':
        with open(file, 'r') as fi:
            while True:
                content = fi.readline()
                if not content:
                    break
                if c := content.rstrip('\n'):
                    seq_list.append(c.split(','))
    elif fmt == 'xlsx':
        import pandas as pd
        df = pd.read_excel(file)
        for i in range(len(df)):
            seq_list.append([*df.iloc[i, :3]])
    else:
        print(f'format {fmt} not recognized!')
        exit()
    return seq_list


def calc_single_content(seq, n1, n2):
    n1_count = seq.count(n1)
    n2_count = seq.count(n2)
    content = (n1_count + n2_count) / len(seq)
    return content


def calc_content(input_path, seq_path, n1='G', n2='C', fmt='t', output_path=None):
    seq_list = read_file(input_path, fmt)
    full_seq = None
    for record in SeqIO.parse(seq_path, 'fasta'):
        full_seq = record.seq.upper()
    if not full_seq:
        print('No sequence captured!')
        return
    if not output_path:
        for name, start, end in seq_list:
            print('\t'.join([name, start, end, f'{calc_single_content(full_seq[int(start) - 1:int(end)], n1, n2):.2f}']))
    else:
        with open(output_path, 'w') as fo:
            for name, start, end in seq_list:
                fo.write('\t'.join(
                    [name, start, end, f'{calc_single_content(full_seq[int(start) - 1:int(end)], n1, n2):.2f}']) + '\n')


if __name__ == '__main__':
    print(sys.argv, file=sys.stderr)
    calc_content(*sys.argv[1:])
