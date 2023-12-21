"""
This script will autocheck every subcommand in ChIP_seq module using test files.
It is generally a normal analysis workflow of ChIP-seq data
Note: The sample data provided are not real ChiP-seq data, and thus could not yield any peaks.
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
    cmd1 = f'wanglab cutadapt -d . -t chip-seq'
    test('cutadapt', cmd1)

    # 2. bowtie
    cmd2 = f'wanglab bowtie -d . -r ./ref.fa -@ 4 -t chip-seq'
    test('bowtie', cmd2)

    # 3. macs3 callpeak
    cmd3 = f'wanglab macs -t ./treat.bam -c control.bam -l eib202 --keep-dup all'
    test('macs3 callpeak', cmd3)


if __name__ == '__main__':
    main()
