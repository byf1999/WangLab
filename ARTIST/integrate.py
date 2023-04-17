'''
This script integrates codes running pipeline of EL-ARTIST and Con-ARTIST.
The functions used are written in the Tn5_function.py script. Ensure they are both in a directory named ARTIST.

Usage and dependencies refers to the EL-ARTIST and Con-ARTIST script.

v1.0 2023.02.21 (byf)
'''
from ARTIST.Tn5_function import *
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
# sam_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\A-V-Tn55_combined_R1_ART.SAM'
window_size = 100
# all_Tn5_sum, Tn5_NT, Tn5_annotations = SAMreader_Tn5(sam_file, chr_search_term, genome_length, chr_start, chr_end, chr_name, window_size=window_size)

## if you have run the SAMreader_Tn5.py script and want to load from .wig2 file, run the following lines
## otherwise, skip it
wig2_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\A-V-Tn55_combined_R1_ART_window_size=100.wig2'
all_Tn5_sum, Tn5_NT, Tn5_annotations = SAMreader_Tn5(wig2_file, chr_search_term, genome_length, chr_start, chr_end,
                                                     chr_name, window_size=window_size)

## now check your data under the specified window_size
saturation = 1 - all_Tn5_sum.count(0) / len(all_Tn5_sum)
print(f'{saturation * 100:.2f}%')

# 3. Sequencing Saturation and Bottleneck Effects


# 4. get unique names
unique_names, unique_indices = get_unique_names(Tn5_annotations)

## check the replication bias
y_max = 25000  # range of # insertions
plt.scatter(Tn5_NT, all_Tn5_sum)
plt.ylim(0, y_max)
plt.show()

# 5. EL-ARTIST
## normalize reads if needed
norm_flag = True # True if normalization is needed, otherwise False
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
cutoffs = calc_cutoffs(all_Tn5_sum)  # recommended cutoffs (25th, 75th percentile and outlier value)
#cutoffs = [2,15,1000]
disc_seq = discretize(all_Tn5_sum, 0, *cutoffs)
## hmm model
essential_calls, count = hmm_essential(disc_seq, regions)
#essential_calls, count = hmm_essential(list(np.squeeze(v)-1), r)
#essential_calls, count = hmm_essential(list(np.squeeze(v)-1), regions)
unique_calls, out_stats = output_in_cds(essential_calls, unique_indices)
print(f'In total {sum(out_stats)} regions:\n'
      f'{out_stats[1]} ({out_stats[1] / sum(out_stats) * 100:.2f}%) are essential\n'
      f'{out_stats[2]} ({out_stats[2] / sum(out_stats) * 100:.2f}%) are domain-essential\n'
      f'{out_stats[0]} ({out_stats[0] / sum(out_stats) * 100:.2f}%) are non-essential')
# sum(essential_calls)

## output to csv
file_path = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\output_test.tsv'
kwargs = {'samreader': [window_size, f'{saturation * 100:.2f}%'],
         'norm': [norm_flag, window_len] if norm_flag else [norm_flag],
         'sw': [sw_size, sim_times, p_value],
         'discritize': [[0, *cutoffs]]}
output_to_tsv(file_path, unique_names, unique_calls, out_stats, **kwargs)

# 6. Con-ARTIST
## read data under experiment condition, the same as above
window_size = 100
wig2_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\A-V-Tn55_combined_R1_ART_window_size=100.wig2'
experiment, *_ = SAMreader_Tn5(wig2_file, chr_search_term, genome_length, chr_start, chr_end,
                                                     chr_name, window_size=window_size)
if norm_flag:
    window_len = 100000
    normed = window_average(Tn5_NT, experiment, window_len, genome_length)
    experiment = list(normed)
## let all_Tn5_sum be the control
control = all_Tn5_sum
## Simulate bottleneck
num_boots = 100
saturation_ratio = experiment.count(0) / control.count(0)
control_norm = simulate_equal_saturation(saturation_ratio, sum(experiment), control, num_boots)
## whether to use HMM
HMM_flag = True
p_value_threshold = 0.05
if HMM_flag:
    fc_thresholds = [-3, -1, 1, 3]
    min_sl, max_sl, min_enr, max_enr, prop_ignore = 0.5, 0.9, 0.5, 0.9, 0.1
    denominators = calc_cutoffs(experiment)
    state_confidence = hmm_train_mwu(experiment, control_norm, essential_calls, unique_indices, p_value_threshold, denominators, fc_thresholds)
    output_cond, call_stats = call_conditional(state_confidence, unique_indices, min_sl, max_sl, min_enr, max_enr, prop_ignore)
else:
    sig_proportion = 0.9
    mwu_boots = run_all_MWU_boots(control_norm, experiment, unique_indices)
    mwu_stats = MWU_summary(mwu_boots, p_value_threshold, sig_proportion, unique_indices, control_norm, experiment)