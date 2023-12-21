import os
import sys
import argparse

__FORMAT__ = {'fasta': ['fasta', 'fa', 'fas', 'fna', 'seq'],
              'xlsx': ['xlsx', 'xls'],
              'txt': ['txt'], }


def merge(args):
    file_path, output_path, all = args.dir, args.out, args.all
    if args.fmt.lower() in __FORMAT__.get('fasta'):
        count = merge_fasta(file_path, output_path, all)
    elif args.fmt.lower() in __FORMAT__.get('xlsx'):
        count = merge_xlsx(file_path, output_path, all)
    elif args.fmt.lower() in __FORMAT__.get('txt'):
        count = merge_txt(file_path, output_path, all)
    else:
        print(f'Format {fmt} not recognized, please check.')
        return

    print(f'Task finished, a total of {count} files are merged.')
    return


def prepare_argparser(argparser):
    argparser.add_argument("-d", "--dir", dest="dir", type=str, required=True, nargs="?",
                           help="Directory of files")
    argparser.add_argument("-f", "--fmt", dest="fmt", type=str, required=False, nargs="?",
                           help="Format of files to be merged, default is fasta", default='fasta')
    argparser.add_argument("-o", "--out", dest="out", type=str, required=False, nargs="?",
                           help="Output file path, default is stdout", default=None)
    argparser.add_argument("-a", "--all", dest="all", type=str, required=False, nargs="?",
                           help="merge all files or check file format in the given path, default is False", default='0')
    return argparser


def merge_fasta(file_path, output_path, all):
    from Bio import SeqIO
    file_list = [file for file in os.listdir(file_path) if os.path.isfile(os.path.join(file_path, file))]
    count = 0
    if output_path:
        fo = open(output_path, 'w')
        for file in file_list:
            if eval(all) or check_file_format(file, 'fasta'):
                count += 1
                print(f'Dealing with {file}...\n', file=sys.stderr)
                for record in SeqIO.parse(os.path.join(file_path, file), 'fasta'):
                    fo.write(f'>{record.id}\n')
                    fo.write(f'{record.seq}\n')
                print(f'Finished writing {file}...\n', file=sys.stderr)
            else:
                print(f'Skip {file}.\n', file=sys.stderr)

        fo.close()
    else:
        for file in file_list:
            if eval(all) or check_file_format(file, 'fasta'):
                count += 1
                print(f'Dealing with {file}...\n', file=sys.stderr)
                for record in SeqIO.parse(os.path.join(file_path, file), 'fasta'):
                    print(record.id)
                    print(record.seq)
                print(f'Finished writing {file}...\n', file=sys.stderr)
            else:
                print(f'Skip {file}.\n', file=sys.stderr)
    return count


def merge_xlsx(file_path, output_path, all):
    import pandas as pd
    file_list = [file for file in os.listdir(file_path) if os.path.isfile(os.path.join(file_path, file))]
    count = 0
    if not output_path:
        print(f'For xlsx format, a output file path must be specified.')
        return
    if eval(all):
        print(f'For xlsx format, every files in the input path will be checked.')
    with pd.ExcelWriter(output_path) as writer:
        for file in file_list:
            if check_file_format(file, 'xlsx'):
                count += 1
                print(f'Dealing with {file}...\n', file=sys.stderr)
                df = pd.read_excel(os.path.join(file_path, file), sheet_name=None)
                for sheet_name in list(df):
                    ddf = pd.read_excel(os.path.join(file_path, file), sheet_name=sheet_name)
                    ddf.to_excel(writer, sheet_name=os.path.splitext(file)[0] + f'_{sheet_name}')
                print(f'Finished writing {file}...\n', file=sys.stderr)
            else:
                print(f'Skip {file}.\n', file=sys.stderr)

    return count


def merge_txt(file_path, output_path, all):
    file_list = [file for file in os.listdir(file_path) if os.path.isfile(os.path.join(file_path, file))]
    count = 0
    if output_path:
        fo = open(output_path, 'w')
        for file in file_list:
            if eval(all) or check_file_format(file, 'txt'):
                count += 1
                print(f'Dealing with {file}...\n', file=sys.stderr)
                with open(os.path.join(file_path, file)) as fi:
                    while True:
                        content = fi.readline()
                        if not content:
                            break
                        fo.write(content)
                print(f'Finished writing {file}...\n', file=sys.stderr)
            else:
                print(f'Skip {file}.\n', file=sys.stderr)

        fo.close()
    else:
        for file in file_list:
            if eval(all) or check_file_format(file, 'txt'):
                count += 1
                print(f'Dealing with {file}...\n', file=sys.stderr)
                with open(os.path.join(file_path, file)) as fi:
                    while True:
                        content = fi.readline()
                        if not content:
                            break
                        print(content)
                print(f'Finished writing {file}...\n', file=sys.stderr)
            else:
                print(f'Skip {file}.\n', file=sys.stderr)

    return count


def check_file_format(file, fmt):
    for f in __FORMAT__.get(fmt):
        if file.endswith(f):
            return True
    return False


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Merge all files in specific format')
    argparser = prepare_argparser(argparser)
    args = argparser.parse_args()
    merge(args)
