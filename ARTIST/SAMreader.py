'''
SAMreader.py
This is a software that can be used to analyze .SAM file.
This method is based on SAMreader_Tn5.m implemented by MatLab in this article:
    Pritchard, J. R., Chao, M. C., Abel, S., Davis, B. M., Baranowski, C., Zhang, Y. J., Rubin, E. J., & Waldor, M. K. (2014). ARTIST: high-resolution genome-wide assessment of fitness using transposon-insertion sequencing. PLoS genetics, 10(11), e1004782. https://doi.org/10.1371/journal.pgen.1004782

This method is implemented by Python and can replace the original one owing to its efficiency (more than 100 times)
The result is a file with .wig2 extension, each line represent the sum of insertions of all sites in the window.
The .wig2 file can be loaded in MatLab as All_Tn5_sum variable to perform further HMM analysis.

Usage:
    Open this script, rectify some parts in the main function(lines under "def main()")
    1. "chr_identifier" is a character specifying the chromosomes in the .sam file,
    usually the same as the term in the first column of tff/gtf file;
    2. "genome_length" is the length of your genome;
    3. "window_size" is the size of window you want to sum reads;
    Finally, move all the .sam file you want to count and this script in a same directory,
        type "python SAMreader.py" or run this script in an IDE such as Pycharm or Spyder.

SAMreader.py v1.0 20230209
v1.1 2023.2.21
    Modify some part of the manual.
'''
import collections
import os
import time


def sam2wig(sam_file, chr_identifier, genome_length, window_size=100):
    all_pos = []
    with open(sam_file, 'r') as infile:
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
    all_pos0 = []
    for i in range(genome_length):
        all_pos0.append(c[i])

    with open(os.path.splitext(sam_file)[0] + f'_window_size={window_size}.wig2', 'w') as outfile:
        for i in range(1, genome_length + 1, window_size):
            outfile.write(f'{sum(all_pos0[i - 1:i - 1 + window_size])}\n')


def main():
    chr_identifier = 'Contig00001'
    genome_length = 4703168
    window_size = 100

    file_list = os.listdir(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM')
    for sam_file in file_list:
        if os.path.splitext(sam_file)[-1] in ('.SAM', '.sam'):
            t0 = time.time()
            print(f'Dealing with {sam_file}...')
            sam2wig(sam_file, chr_identifier=chr_identifier, genome_length=genome_length, window_size=window_size)
            print(f'Finished, cost {time.time() - t0:.2f} seconds\n')


# sam_file = 'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\A-V-Tn55_combined_R1_ART.SAM'
if __name__ == '__main__':
    main()
