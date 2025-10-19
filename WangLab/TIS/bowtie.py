import os
import sys
import time


def bowtie(dir_path, sfx, num_threads=2):
    assert os.path.exists(dir_path), "Can not open this directory";

    file_list = os.listdir(dir_path)
    for f in file_list:
        if f.endswith(sfx):
            t = time.localtime()
            print(f"\n\n===========Beginning at {t.tm_hour:0>2}:{t.tm_min:0>2}:{t.tm_sec:0>2}===============\n")
            file = os.path.join(dir_path, f)
            name = f.replace(sfx, '.sam')
            out = os.path.join(dir_path, name)
            print(f"Fastq file -> {file}\n")
            print(f"Output file-> {out}\n\n")
            cmd = f'bowtie -v 3 -a --best --strata -m 1 -q --sam --threads {num_threads} index {file} {out} 2> {os.path.splitext(out)[0]}.log'
            print(cmd)
            os.system(cmd)
            print(f"\n---------- Finished for: {file} ----------------\n\n")


def main(args):
    l = len(args)
    if l not in (3, 4):
        print(f"Usage: python {args[0]} 1.dir_trimmed files 2.Reference(../EIB202) 3.Num_threads[optional]\n")
        return

    sfx = "_trimmed.fq"
    cmd_build_index = f'bowtie-build -f {args[2]} index'
    # print(cmd_build_index)
    os.system(cmd_build_index)
    bowtie(args[1], sfx, args[3] if l == 4 else 1)


if __name__ == '__main__':
    main(sys.argv)
