import argparse
import sys
import os
import time


def prepare_argparser(argparser):
    argparser.add_argument("-i", "--input", dest="input", type=str, required=True, nargs="?",
                           help="Path of reference sequence")
    argparser.add_argument("-d", "--dir", dest="dir", type=str, required=True, nargs="?",
                           help="Input file path")
    argparser.add_argument("-s", "--sfx", dest="sfx", type=str, required=False, nargs="?",
                           help="Suffix of input files, default is .fq.gz", default='.fq.gz')
    argparser.add_argument("--remove-tmp", dest="remove_tmp", type=str, required=False, nargs="?", default=True,
                           help="Whether to remove temporary files, include trimmed, sam and unsorted bam, default is yes")
    argparser.add_argument("-p", "--threads", dest="threads", type=str, required=False, nargs="?", default='1',
                           help="Number of threads used, default is 1")

    return


def run(args):
    dir = args.dir
    sfx = args.sfx
    remove_temp = args.remove_tmp
    threads = args.threads

    cmd_build_index = f'bowtie-build -f {args.input} index'
    # print(cmd_build_index)
    os.system(cmd_build_index)
    file_list = os.listdir(dir)
    for file in file_list:
        if not sfx or file.endswith(sfx):
            file = cutadapt(os.path.join(dir, file), sfx)
            bowtie(file, remove_temp, threads)
            # htseq2(file)


def cutadapt(file, sfx):
    # adp = "GATCGGAAGAGCACACGTCT"
    adp = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    pm = " -O 5 -e 0.1 -m 8"
    print(f"\n\n================Trimming {sfx} adapter ====================\n")
    out = file.replace(sfx, '_trimmed.fq')
    print(f"Fastq file ->{file}\n")
    print(f"output file->{out}\n\n")
    cmd = f"cutadapt -a {adp} {pm} --cores 0 -o {out} {file} > {os.path.basename(file).split('.')[0]}.log"
    print(cmd)
    os.system(cmd)
    print(f"\n---------- Trimming ended for: {file} ----------------\n\n")
    return out


def bowtie(file, remove_temp, threads):
    t = time.localtime()
    print(f"\n\n===========Beginning at {t.tm_hour:0>2}:{t.tm_min:0>2}:{t.tm_sec:0>2}===============\n")
    out = file.replace('_trimmed.fq', '.sam')
    print(f"Fastq file -> {file}\n")
    print(f"Output file-> {out}\n\n")
    cmd = f"bowtie -M 1 -q --sam -p {threads} index {file} {out} 2>> {os.path.basename(file.replace('_trimmed.fq', '')).split('.')[0]}.log"
    print(cmd)
    os.system(cmd)
    os.system(f"samtools view -bSh {out} -@ {int(threads) - 1} > {out.replace('.sam', '_tmp.bam')}")
    os.system(f"samtools sort {out.replace('.sam', '_tmp.bam')} -@ {int(threads) - 1} > {out.replace('sam', 'bam')}")
    os.system(f"samtools index {out.replace('sam', 'bam')} -@ {threads}")
    print(f"\n---------- Finished for: {file}----------------\n\n")
    if remove_temp:
        if sys.platform in ('linux', 'darwin'):
            rm_cmd = 'rm'
        elif sys.platform == 'win32':
            rm_cmd = 'del'
        else:
            return
        os.system(f"{rm_cmd} {file} {out} {out.replace('.sam', '_tmp.bam')}")


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    prepare_argparser(argparser)
    args = argparser.parse_args()
    run(args)
