'''
This script include workflow of EL-ARTIST pipeline dealing with Tn5-based transposon.
The functions used are written in the Tn5_function.py script. Ensure they are both in a directory named ARTIST.

Usage: run the code line by line based on the notes.

Dependency: matplotlib, numpy, pandas, scikit-learn

v1.0 2023.02,21 (byf)
'''
from ARTIST.Tn5_function import *
from ARTIST.check_seq_depth import check_seq_depth
import matplotlib.pyplot as plt

# 1. genome_parser
gtf_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\genomic_gene.gff'
chr_search_term = 'Contig00001'
genome_length = 4703168

cds_search_term = r'ID=(.+?);'
## cds_search_term choices:
## r'ID=(.+?);' for ID=GE000001;
## r'gene_id "(.+?)"' for gene_ID "GE000001";
chr_start, chr_end, chr_name = genome_parser(gtf_file, chr_search_term, genome_length, cds_search_term=cds_search_term)

# 2. SAMreasder
window_size = 30
input_file = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\wig2\WT_window_size=100.wig2'
all_Tn5_sum, Tn5_NT, Tn5_annotations = SAMreader_Tn5(input_file, chr_search_term, genome_length, chr_start, chr_end,
                                                     chr_name, window_size=window_size)

## now check your data under the specified window_size
ax = check_seq_depth(all_Tn5_sum, sep=1e5)

ax.get_figure().savefig(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM\pics\depth.svg')
saturation = 1 - all_Tn5_sum.count(0) / len(all_Tn5_sum)
print(f'{saturation * 100:.2f}%')

# 3. Sequencing Saturation and Bottleneck Effects


# 4. get unique names
unique_names, unique_indices = get_unique_names(Tn5_annotations)

## check the replication bias
y_max = 25000  # range of # insertions
plt.scatter(Tn5_NT, all_Tn5_sum)
plt.ylim(0, y_max)
plt.xlabel('# of reads')
plt.ylabel('# of unique insertions')
plt.show()

# 5. EL-ARTIST
## normalize reads if needed
norm_flag = True  # True if normalization is needed, otherwise False
if norm_flag:
    window_len = 100000
    normed = window_average(Tn5_NT, all_Tn5_sum, window_len, genome_length)
    all_Tn5_sum = list(normed)
## Sliding window
p_value = 0.08
sim_times = 1000
sw_size = 7
regions, values = sliding_window(all_Tn5_sum, p_value, sim_times, sw_size)
print(f'the percent of initial essential is {float(sum(regions) / len(values) * 100):.2f}%')
### 0 for non-essential and 1 for essential
### in python ,there is no need to add 1 to regions

### rectify p_value based on values and re-run sliding_window

## discretize
cutoffs = calc_cutoffs(all_Tn5_sum)
cutoffs = [1, 1, 5]
### recommended cutoffs (25th, 75th percentile and outlier value)
### you can also specify cutoffs as you like
### eg. cutoffs = [1, 5, 17]
disc_seq = discretize(all_Tn5_sum, 0, *cutoffs)
## hmm model
essential_calls, count = hmm_essential(disc_seq, regions)
unique_calls, out_stats = output_in_cds(essential_calls, unique_indices)
print(f'In total {sum(out_stats)} regions:\n'
      f'{out_stats[1]} ({out_stats[1] / sum(out_stats) * 100:.2f}%) are essential\n'
      f'{out_stats[2]} ({out_stats[2] / sum(out_stats) * 100:.2f}%) are domain-essential\n'
      f'{out_stats[0]} ({out_stats[0] / sum(out_stats) * 100:.2f}%) are non-essential')

## output to csv
file_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\ARTIST\EL-ARTIST\output_2.tsv'
kwargs = {'fi_path': [input_file],
          'samreader': [window_size, f'{saturation * 100:.2f}%'],
          'norm': [norm_flag, window_len] if norm_flag else [norm_flag],
          'sw': [sw_size, sim_times, p_value],
          'discritize': [[0, *cutoffs]]}
output_to_tsv(file_path, unique_names, unique_calls, out_stats, method='EL-ARTIST', **kwargs)

output_essential_calls = True
file_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\ARTIST\EL-ARTIST\final_essential_2.tsv'
if output_essential_calls:
    output_to_tsv(file_path, essential_calls, essential_calls, out_stats, method='EL-ARTIST', **kwargs)