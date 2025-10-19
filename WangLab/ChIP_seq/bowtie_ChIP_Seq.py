"""
This script is for ChIP-seq.
"""
import os
import sys
import time


def bowtie(dir_path, sfx, num_cores=2):
    assert os.path.exists(dir_path), "Can not open this directory";

    file_list = os.listdir(dir_path)
    for f in file_list:
        if f.endswith(sfx):  # for ChIP-seq, read1 only
            t = time.localtime()
            print(f"\n\n===========Beginning at {t.tm_hour:0>2}:{t.tm_min:0>2}:{t.tm_sec:0>2}===============\n")
            file = os.path.join(dir_path, f)
            out = f.replace(sfx, '.sam')
            print(f"Fastq file -> {file}\n")
            print(f"Output file-> {out}\n\n")
            cmd = f'bowtie -M 1 -q --sam -p {num_cores} index {file} {out} 2>  {os.path.splitext(out)[0]}.log'
            print(cmd)
            os.system(cmd)
            os.system(f"samtools view -bSh {out} -@ {int(num_cores) - 1} > {out.replace('.sam', '_tmp.bam')}")
            os.system(f"samtools sort {out.replace('.sam', '_tmp.bam')} -@ {int(num_cores) - 1} > {out.replace('sam', 'bam')}")

            os.system(f"samtools index {out.replace('sam', 'bam')} -@ {num_cores}")
            print(f"\n---------- Finished for: {file}----------------\n\n")


def main(args):
    l = len(args)
    if not l in (3, 4):
        print(f"Usage: python {args[0]} 1.dir_trimmed files 2.Reference(../EIB202) 3.Num_cores[optional]\n")
        return

    sfx = "_trimmed.fq"
    cmd_build_index = f'bowtie-build -f {args[2]} index'
    # print(cmd_build_index)
    os.system(cmd_build_index)
    bowtie(args[1], sfx, args[3] if l == 4 else 2)

    if sys.platform in ('linux', 'darwin'):
        rm_cmd = 'rm'
    elif sys.platform == 'win32':
        rm_cmd = 'del'
    else:
        return
    os.system(f'{rm_cmd} {os.path.join(args[1], ("*" + "_tmp.bam"))}')


if __name__ == '__main__':
    main(sys.args)
