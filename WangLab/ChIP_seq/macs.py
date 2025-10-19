"""
This script contain several default parameters setting of MACS3.
Multiple files can be entered at the same time in the format of <file1,file2,file3...>
Each treatment file will be analyzed with every control files in turn.

Usage: python macs.py <treatment_files> <control_files> <bacteria_name/genome_length>

version 1.0 20230608 by Byf
"""
import os
import sys


def main(args, keepduplicates='1'):
    if not len(args) in (4, 5):
        print('Usage: python macs.py <treatment_files> <control_files> <bacteria_name/genome_length>')
        return

    genome_length_dict = {'eib202': 3760363,
                          'A_veronii': 4703168, }

    test_files = args[1].split(',')  # treatment file
    control_files = args[2].split(',')  # control file
    genome_length = int(genome_length_dict.get(args[3], args[3]))  # genome length
    q_value = 0.01  # q_value threshold
    MFold = [1, 50]

    for test in test_files:
        t = os.path.basename(test)
        for control in control_files:
            c = os.path.basename(control)
            out = f'{t.split("_")[0]}{c.split("_")[0]}'
            cmd = f'macs3 callpeak -t {test} -c {control} -f BAM --outdir {out} -n {out} -g {genome_length} ' \
                  f'-q {q_value} -B -m {MFold[0]} {MFold[1]} --keep-dup {keepduplicates} 2> {out}.log'
            print(cmd)
            os.system(cmd)


if __name__ == '__main__':
    main(sys.argv)
