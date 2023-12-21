"""
This script will autocheck every subcommand in RNA_seq module using test files.
It is generally a normal analysis workflow of RNA-seq data
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
    os.chdir('./input_files')
    # 1. cutadapt
    cmd1 = f'wanglab cutadapt -d . -p Tn5 -t tn-seq'
    test('cutadapt', cmd1)

    # 2. bowtie
    cmd2 = f'wanglab bowtie -d . -r ./ref.fa -@ 4 -t tn-seq'
    test('bowtie', cmd2)

    # 3. count_reads
    cmd3 = f'wanglab count_reads -i ./annot.gff -c Contig00001 -r "ID=(.+?);" -l 4703168 -d .'
    test('count_reads', cmd3)

    # 4. combine_reads
    cmd4 = f'wanglab combine_reads -d .'
    test('combine_reads', cmd4)


if __name__ == '__main__':
    main()
