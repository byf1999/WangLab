"""
This script is for ChIP-seq data.
"""
##version:
# v0.1 add "rm tmp"
# v0.41 add 2 base to 5' adapter, set -O 5 -e 0.1. proved to be better than v0.3
# v0.5 change the 5' adapter to 19 bp
# v0.54 solve 2 bug because of cutadapt update -by byf
# v1.0 implemented by Python, remove unzip, directly use .gz file
#      auto-detect cpu cores to be used.
import os
import sys


def cutadapt(ip, s, adp, op, pwd, pm):
    file_list = os.listdir(pwd)
    for f in file_list:
        if f.endswith(ip):
            print(f"\n\n================Trimming {ip} adapter ====================\n")
            file = os.path.join(pwd, f)
            out = file.replace(ip, op)
            print(f"Fastq file ->{file}\n")
            print(f"output file->{out}\n\n")
            cmd = f"cutadapt -{s} {adp} {pm} --cores 0 -o {out} {file} >> cutadapt.log"
            print(cmd)
            os.system(cmd)
            print(f"\n---------- Trimming ended for: {file} ----------------\n\n")


def main(args):
    if not len(args) == 2:
        print(f"Usage: python {arggs[0]} <dir_of_fq.gz_file>")
        return
    # for ChIP-seq, 3' trimming only.
    ip = '.fq.gz'  # input file suffix
    str = "a"
    adp = "GATCGGAAGAGCACACGTCT"
    op = "_trimmed.fq"
    pm = " -O 5 -e 0.1 -m 8"
    pwd = args[1]
    cutadapt(ip, str, adp, op, pwd, pm)
    print("\n****************trimming finished************************\n")


if __name__ == '__main__':
    main(sys.argv)
