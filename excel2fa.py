"""
This script transforms a xlsx/xls file to fasta format.
The xlsx/xls file contains two columns, the first represents gene names and the second is the corresponding sequences.


Usage: python excel2fa.py <excel_file> [output_file]
"""

import os
import pandas as pd
import sys


def exception(*args, **kwargs):
    raise Exception('The format of Input file is not supported!')


def main():
    assert len(sys.argv) in (2, 3), "Usage: python excel2fa.py <excel_file> [output_file]"

    fi_path = sys.argv[1]
    fi, ext = os.path.splitext(fi_path)
    fo_path = sys.argv[2] if len(sys.argv) == 3 else fi + '.fa'
    if os.path.exists(fo_path):
        print(f'The file {fo_path} already exists!')
        return

    fn_dict = {'.xlsx': pd.read_excel, '.xls': pd.read_excel, '.csv': pd.read_csv, '.tsv': pd.read_table}

    df = fn_dict.get(ext, exception)(fi_path, header=0)  # header line exists.
    with open(fo_path, 'w') as fo:
        for i in range(len(df)):
            fo.write('>' + df.iloc[i, 0] + '\n')
            fo.write(df.iloc[i, 1] + '\n')


if __name__ == '__main__':
    main()
