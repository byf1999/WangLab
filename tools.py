import collections
import os


def sam2wig(sam_file, chr_identifier, file_type='sam'):
    if file_type == 'sam':
        file_type = 'r'
    elif file_type == 'bam':
        file_type = 'rb'
    all_pos = []
    with open(sam_file, file_type) as infile:
        while True:
            content = infile.readline()
            if not content.startswith('@'): break
        while content:
            content = content.strip('\n').split('\t')
            if not content[2] == chr_identifier:
                content = infile.readline()
                continue
            if content[1] == '0':
                pos = int(content[3])
            elif content[1] == '16':
                pos = int(content[3]) + len(content[9]) - 2
            else:
                content = infile.readline()
                continue
            all_pos.append(pos)
            content = infile.readline()

    c = collections.Counter(all_pos)

    with open(os.path.splitext(sam_file)[0] + '.wig', 'w') as outfile:
        for k in sorted(c): # only inserted sites
            outfile.write(f'{k}\t{c[k]}\n')
