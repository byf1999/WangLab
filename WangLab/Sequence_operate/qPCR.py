import sys
import os
import argparse

import numpy as np
import pandas as pd


def prepare_argparser(argparser):
    argparser.add_argument("-i", '--input', dest='input', type=str, required=True, nargs="+",
                           help='Input file, the original xls file generated by QS3 software')
    argparser.add_argument("-s", '--sample', dest='ctrl_sample', type=str, required=True, nargs="?",
                           help='Control sample name, should be the same in xls file')
    argparser.add_argument("-t", '--target', dest='ctrl_target', type=str, required=True, nargs="?",
                           help='Reference gene name, should be the same in xls file')
    argparser.add_argument("-o", '--output', dest='output', type=str, required=False, nargs="?",
                           help='Output file path, if not specified, a xlsx file with same name will be generated.')
    argparser.add_argument("-f", '--force', dest='force', type=str, required=False, nargs='?',
                           help='Set -f True means force to overwrite output file when it\'s already existed. Default is False',
                           default=False)


def calc(args):
    input_files = args.input
    if isinstance(input_files, list):
        for input_file in input_files:
            args.input = input_file
            calc(args)
    elif isinstance(input_files, str):
        if os.path.isfile(input_files):
            args.input = input_files
            calc_single_file(args)
        else:
            for input_file in os.listdir(input_files):
                args.input = input_file
                calc_single_file(args)
    else:
        print('No input file caught.')


def try_calc(func):
    def wrapper(args):
        input_file = args.input
        print(f'Dealing with {input_file}...')
        # try:
        func(args)
        print('Finished!\n')
        # except Exception as e:
        #     print(f'Failed during {input_file}, with exception error:\n{e}\nContinue to solve the next file.\n')

    return wrapper


@try_calc
def calc_single_file(args):
    """
    # This is test data
    input_file = '/home/byf1999/Documents/Tn-seq/qPCR/2023-09-14_150158.xls'
    ctrl_sample = '-1-1'
    ctrl_target = 'gyrB'
    """
    ctrl_sample = args.ctrl_sample
    ctrl_target = args.ctrl_target
    output_file = args.output
    input_file = args.input

    df = pd.read_excel(input_file, sheet_name='Results', header=43, dtype=str)

    # Check if Sample name and target name is correct.
    if not df['Sample Name'].isin([ctrl_sample]).any():
        print(f"{ctrl_sample} not found in sample names, please check!")
        return
    if not df['Target Name'].isin([ctrl_target]).any():
        print(f"{ctrl_target} not found in target names, please check!")
        return

    data_dict = {}
    for row, (sample, target, ct) in enumerate(zip(df['Sample Name'], df['Target Name'], df['CT'])):
        # Deal with Undetermined values
        v = float(ct) if not ct == 'Undetermined' else np.nan
        i = 1
        while (f'{sample}_{i}', target) in data_dict:
            # Deal with biological replicate
            i += 1
        data_dict[(f'{sample}_{i}', target)] = [v]
        df['Sample Name'][row] = f'{sample}_{i}'

    for (sample, _), value in data_dict.items():
        value.append(value[0] - data_dict[(sample, ctrl_target)][0])  # d_ct = ct - ct_reference_gene

    for (_, target), value in data_dict.items():
        value.append(value[1] - data_dict[(f'{ctrl_sample}_1', target)][1])  # dd_Ct = d_ct - d_ct_control

    d_ct = []
    dd_ct = []
    _2_dd_ct = []
    for sample, target, ct in zip(df['Sample Name'], df['Target Name'], df['CT']):
        value = data_dict[(sample, target)]
        d_ct.append(value[1])
        dd_ct.append(value[2])
        _2_dd_ct.append(2 ** (-value[2]))

    df['d_Ct'] = d_ct
    df['dd_Ct'] = dd_ct
    df['2^-dd_Ct'] = _2_dd_ct

    output_file = output_file if output_file else input_file + 'x'
    if (not args.force) and os.path.exists(output_file):
        print(f'{output_file} already exists!')
    else:
        dfv = []
        samples = sorted(list(set(df['Sample Name'])), key=list(df['Sample Name']).index)
        targets = sorted(list(set(df['Target Name'])), key=list(df['Target Name']).index)
        for sample in samples:
            dfv.append([])
            for target in targets:
                try:
                    tmp = float(df[(df['Target Name'] == target) & (df['Sample Name'] == sample)]['2^-dd_Ct'])
                    dfv[-1].append(tmp)
                except:
                    dfv[-1].append(np.nan)
        pd.DataFrame(dfv, index=samples, columns=targets).to_excel(output_file)

        # else:
        #     df[['Sample Name', 'Target Name', '2^-dd_Ct']].to_excel(output_file, index=False)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    prepare_argparser(argparser)
    args = argparser.parse_args()
    calc(args)