import os
import sys
import time

def bowtie(dir_path, sfx, num_threads=2):
    assert os.path.exists(dir_path), "Can not open this directory";

    file_list = os.listdir(dir_path)
    for f in file_list:
        if f.endswith(sfx):
            t = time.localtime()
            print(f"\n\n===========Beginning at {t.tm_hour}:{t.tm_min}:{t.tm_sec}===============\n")
            file = os.path.join(dir_path, f)
            name = f.replace(sfx, '.sam')
            out = os.path.join(dir_path, name)
            print(f"Fastq file -> {file}\n")
            print(f"Output file-> {out}\n\n")
            cmd = f'bowtie -v 3 -a --best --strata -m 1 -q --threads {num_threads} index {file} {out}'
            print(cmd)
            os.system(cmd)
            print(f"\n---------- Finished for: {file} ----------------\n\n")


l = len(sys.argv)
assert l in (3, 4), f"Usage: python {sys.argv[0]} 1.dir_trimmed files 2.Reference(../EIB202)\n"

sfx = "_trimmed.fq"
cmd_build_index = f'bowtie-build -f {sys.argv[2]} index'
# print(cmd_build_index)
os.system(cmd_build_index)
bowtie(sys.argv[1], sfx, sys.argv[3] if l == 4 else 2)
