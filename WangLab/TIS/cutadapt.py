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


def cutadapt(ip, s, adp, op, pwd, pm):
    file_list = os.listdir(pwd)
    for f in file_list:
        if f.endswith(ip):
            t = time.localtime()
            print(f"\n\n================Trimming {ip} adapter ====================\n")
            print(f"\n\n===========Beginning at {t.tm_hour}:{t.tm_min}:{t.tm_sec}===============\n")
            file = os.path.join(pwd, f)
            name = f.replace(ip, '') + op
            out = os.path.join(pwd, name)
            print(f"Fastq file ->{file}\n")
            print(f"output file->{out}\n\n")
            cmd = f"cutadapt -{s} {adp} {pm} --cores 0 -o {out} {file} >> cutadapt.log"
            print(cmd)
            os.system(cmd)
            print(f"\n---------- Trimming ended for: {file} ----------------\n\n")


adp_dict = {'TN5': "AAGCTTCGGCCGCCTAGGCC",
            'PSC189': "GACTTATCAGCCAACCTGT",
            'PMAR2XT7': "GACTTATCAGCCAACCTGT", }


def main(args, plasmid='TN5'):
    if not len(args) == 2:
        print(f"Usage: python {sys.argv[0]} <dir_of_fq.gz_file> [plasmid name]\n"
              f"plasmid name include TN5, PSC189, PMAR2XT7")
        return

    # step 1: 5' trimming
    ip = ".fq.gz"
    str = "g"
    # adp = "GACTTATCAGCCAACCTGT" # pSC189
    # adp = "GACTTATCAGCCAACCTGTTA" # pMar2xT7
    # adp = "AAGCTTCGGCCGCCTAGGCC" # pUTmini-Tn5Km2
    if plasmid in adp_dict:
        adp = adp_dict[plasmid.upper()]
    else:
        print(f'Plasmid {plasmid} not recognized!')
        return
    op = "_5_tmp.fq"
    pm1 = "-O 17 -e 0.2  --match-read-wildcards --discard-untrimmed"
    pwd = args[1]
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

    if sys.platform in ('linux', 'darwin'):
        rm_cmd = 'rm'
    elif sys.platform == 'win32':
        rm_cmd = 'del'
    else:
        return
    os.system(f'{rm_cmd} {os.path.join(pwd, ("*" + ip))}')


if __name__ == '__main__':
    main(sys.argv)
