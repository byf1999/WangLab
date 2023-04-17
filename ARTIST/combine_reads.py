import sys
import pandas as pd
import os

assert len(sys.argv) == 2, 'Usage: python combine_reads.py <_reads.csv dir>'
fi_path = sys.argv[1]
df_dict = {}
sfx = '_reads.csv'

file_list = os.listdir(fi_path)
for file in file_list:
    if file.endswith(sfx):
        print(f"Reading {file}...")
        df = pd.read_table(os.path.join(fi_path, file), sep=',', header=None)
        if df.iloc[1, 0] in df_dict:
            df_dict[df.iloc[1, 0]][file.replace(sfx, '')] = df[3]
        else:
            df.columns = ['Feature', 'Start', 'End', file.replace(sfx, '')]
            df_dict[df.iloc[1, 0]] = df

print("Writing...")
for k, df in df_dict.items():
    df.to_excel(os.path.join(fi_path, f'Combination_reads_{k}.xlsx'), index=False)
