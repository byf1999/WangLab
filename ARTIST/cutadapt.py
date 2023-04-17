##version:
# v0.1 add "rm tmp"
# v0.41 add 2 base to 5' adapter, set -O 5 -e 0.1. proved to be better than v0.3
# v0.5 change the 5' adapter to 19 bp
# v0.54 solve 2 bug because of cutadapt update -by byf
# v1.0 implemented by Python, remove unzip, directly use .gz file
#      auto-detect cpu cores to be used.
import os
import sys
import time

assert len(sys.argv) == 2, "Usage: python cutadapt.py <dir_of_fq.gz_file>"


def cutadapt(ip, s, adp, op, pwd, pm):
    file_list = os.listdir(pwd)
    for f in file_list:
        if f.endswith(ip):
            t = time.localtime()
            print(f"\n\n==============s==Trimming {ip} adapter ====================\n")
            print(f"\n\n===========Beginning at {t.tm_hour}:{t.tm_min}:{t.tm_sec}===============\n")
            file = os.path.join(pwd, f)
            name = f.replace(ip, '') + op
            out = os.path.join(pwd, name)
            print(f"Fastq file ->{file}\n")
            print(f"output file->{out}\n\n")
            print(f"cutadapt -{s} {adp} {pm} --cores 0 -o {out} {file}")
            os.system(f"cutadapt -{s} {adp} {pm} --cores 0 -o {out} {file}")
            print(f"\n---------- Trimming ended for: {file} ----------------\n\n")


def main():
    # step 1: 5' trimming
    ip = ".fq.gz"
    str = "g"
    adp = "GACTTATCAGCCAACCTGT"
    op = "_5_tmp.fq"
    pm1 = "-O 17 -e 0.2  --match-read-wildcards --discard-untrimmed"
    pwd = sys.argv[1]
    cutadapt(ip, str, adp, op, pwd, pm1)
    print("\n****************step1 finished************************\n")

    # step 2: 3' trimming
    ip = op
    str = "a"
    adp = "ATACCACGAC"
    op = "_trimmed.fq"
    pm = " -O 5 -e 0.1 -m 10"
    cutadapt(ip, str, adp, op, pwd, pm)
    print("\n****************step2 finished************************\n")

    if sys.platform == 'linux':
        rm_cmd = 'rm'
    elif sys.platform.startswith('win'):
        rm_cmd = 'del'
    else:
        return
    os.system(f'{rm_cmd} {os.path.join(pwd, ("*" + ip))}')


if __name__ == '__main__':
    main()
