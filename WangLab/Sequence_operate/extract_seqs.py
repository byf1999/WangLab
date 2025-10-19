import sys
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord


def read_file(file, fmt):
    seq_list = []
    if fmt == 't':
        with open(file, 'r') as fi:
            while True:
                content = fi.readline()
                if not content:
                    break
                if c := content.rstrip('\n'):
                    flag = len(c.split('\t'))
                    if flag == 3:
                        seq_list.append(c.split('\t'))
                    elif flag == 1:
                        seq_list.append(c)
    elif fmt == ',':
        with open(file, 'r') as fi:
            while True:
                content = fi.readline()
                if not content:
                    break
                if c := content.rstrip('\n'):
                    flag = len(c.split(','))
                    if flag == 3:
                        seq_list.append(c.split(','))
                    elif flag == 1:
                        seq_list.append(c)
    elif fmt == 'xlsx':
        import pandas as pd
        df = pd.read_excel(file)
        flag = 3
        for i in range(len(df)):
            seq_list.append([*df.iloc[i, :flag]])
    else:
        print(f'format {fmt} not recognized!')
        exit()
    return seq_list, flag


__FLAG_DICT__ = {
    "translate": Seq.translate,
    "complement": Seq.complement,
    "complement_rna": Seq.complement_rna,
    "reverse_complement": Seq.reverse_complement,
    "reverse_complement_rna": Seq.reverse_complement_rna,
    "rna": lambda x: x.upper().replace('T', 'U')
}


def extract_content(input_path, seq_path, fmt='t', output_path=None, translate_flag=False):
    def generate_seq(seq, s=0, e=0):
        s, e = int(s), int(e)
        flag_dict = __FLAG_DICT__
        if not (s and e):
            return flag_dict[translate_flag](seq) if translate_flag else seq
        if s <= e:
            return flag_dict[translate_flag](seq[int(s) - 1:int(e)]) if translate_flag else seq[int(s) - 1:int(e)]
        else:
            s, e = e, s
            tmp_seq = Seq.reverse_complement(seq[int(s) - 1:int(e)])
            return flag_dict[translate_flag](tmp_seq) if translate_flag else tmp_seq

    seq_list, flag = read_file(input_path, fmt)
    if flag == 3:
        record = list(SeqIO.parse(seq_path, 'fasta'))[0]
        full_seq = record.seq.upper()
        if not full_seq:
            print('No sequence captured!')
            return
        if not output_path:
            for name, start, end in seq_list:
                print(f'>{name}\n{generate_seq(full_seq, start, end)}')
        else:
            records = []
            for name, start, end in seq_list:
                tmp_seq = generate_seq(full_seq, start, end)
                records.append(SeqRecord(Seq.Seq(tmp_seq, ), id=name, description=''))
            SeqIO.write(records, output_path, 'fasta')
    elif flag == 1:
        if not output_path:
            for record in SeqIO.parse(seq_path, 'fasta'):
                if record.id in seq_list:
                    print(f'>{record.id}\n{record.seq}')
        else:
            records = []
            for record in SeqIO.parse(seq_path, 'fasta'):
                if record.id in seq_list:
                    records.append(
                        SeqRecord(Seq.Seq(generate_seq(record.seq), ), id=record.id, description=record.description))
            SeqIO.write(records, output_path, 'fasta')


if __name__ == '__main__':
    print(sys.argv, file=sys.stderr)
    extract_content(*sys.argv[1:])
