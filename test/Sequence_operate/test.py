"""
This script will autocheck every subcommand in Sequence_operate module using test files.
"""

import os


def test(name, cmd):
    print(f'Testing {name}...')
    try:
        os.system(cmd)
        print('Successed.')
    except Exception:
        print('Failed.')


def main():
    try:
        os.mkdir('tmp')
    except FileExistsError:
        pass
    # Test for extract_seqs
    # extract sequences of given positions
    cmd1 = f'wanglab extract_seqs -i ./input_files/pos.txt -r ./input_files/ref.fa'
    # extract sequence and output to given file
    cmd2 = f'wanglab extract_seqs -i ./input_files/pos.txt -r ./input_files/ref.fa -o ./tmp/out2.fa'
    # extract sequence and output the complement sequence to given file
    cmd3 = f'wanglab extract_seqs -i ./input_files/pos.txt -r ./input_files/ref.fa -o ./tmp/out3.fa -t complement'
    # the position file in comma-separated
    cmd4 = f'wanglab extract_seqs -i ./input_files/pos.csv -r ./input_files/ref.fa -o ./tmp/out4.fa -t complement -f ,'
    [test('extract_seqs', cmd) for cmd in [cmd1, cmd2, cmd3, cmd4]]
    # extract sequences of given seq names
    cmd5 = f'wanglab extract_seqs -i ./input_files/names.txt -r ./input_files/seqs.fa -o ./tmp/out5.fa'
    test('extract_seqs', cmd5)

    # Test for calc_content
    # calculate GC content of given positions
    cmd6 = f'wanglab calc_content -i ./input_files/pos.txt -r ./input_files/ref.fa -n G C'
    # calculate and output to given file
    cmd7 = f'wanglab calc_content -i ./input_files/pos.txt -r ./input_files/ref.fa -n G C -o ./tmp/out7.fa'
    # calculate AT content of a comma-separated input file
    cmd8 = f'wanglab calc_content -i ./input_files/pos.csv -r ./input_files/ref.fa -n A T -o ./tmp/out8.fa -f ,'
    [test('calc_content', cmd) for cmd in [cmd6, cmd7, cmd8]]

    # Test for file_merge
    # merge all files in input_files that ends with 'txt'
    cmd9 = f'wanglab file_merge -d ./input_files -f txt -o ./tmp/out9.txt'
    # merge all files in input_files that ends with 'fa'
    cmd10 = f'wanglab file_merge -d ./tmp -f fa -o ./tmp/out10.txt'
    # merge all files in input_files that ends with 'xlsx
    cmd11 = f'wanglab file_merge -d ./input_files -f xlsx -o ./tmp/out11.xlsx'
    [test('file_merge', cmd) for cmd in [cmd9, cmd10, cmd11]]


if __name__ == '__main__':
    main()
