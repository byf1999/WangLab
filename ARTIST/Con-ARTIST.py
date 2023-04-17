'''
This script include workflow of EL-ARTIST pipeline dealing with Tn5-based transposon.
The functions used are written in the Tn5_function.py script. Ensure they are both in a directory named ARTIST.

Usage: run the code line by line based on the notes.

Dependency: matplotlib, numpy, pandas, scikit-learn

v1.0 2023.02,21 (byf)
'''
from ARTIST.Tn5_function import *

# 1. EL-ARTIST part, to get final calls of the control library
#   Parameters are the same in the EL-ARTIST script
gtf_file = r'F:\LearningFiles\Master\8.TIS\ARTIST\A. veronii 2021\data\genomic_gene.gff'
chr_search_term = 'GD2019'
genome_length = 4703168
cds_search_term = r'ID=(.+?);'
## r'gene_id "(.+?)"'
window_size = 100
input_file = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\wig2\WT_window_size=100.wig2'

chr_start, chr_end, chr_name = genome_parser(gtf_file, chr_search_term, genome_length, cds_search_term=cds_search_term)
all_Tn5_sum, Tn5_NT, Tn5_annotations = SAMreader_Tn5(input_file, chr_search_term, genome_length, chr_start, chr_end,
                                                     chr_name, window_size=window_size)
unique_names, unique_indices = get_unique_names(Tn5_annotations)

## norm
norm_flag = True  # True if normalization is needed, otherwise False
if norm_flag:
    window_len = 100000
    normed = window_average(Tn5_NT, all_Tn5_sum, window_len, genome_length)
    all_Tn5_sum = list(normed)

## EL-ARTIST or just load result from a file.
control_call_file = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\ARTIST\EL-ARTIST\final_essential_1.tsv'
# control_call_file = None
if control_call_file:
    essential_calls = load_essential_calls(control_call_file)
else:
    ## Sliding window
    p_value = 0.03
    sim_times = 1000
    sw_size = 7
    regions, values = sliding_window(all_Tn5_sum, p_value, sim_times, sw_size)
    cutoffs = calc_cutoffs(all_Tn5_sum)
    disc_seq = discretize(all_Tn5_sum, 0, *cutoffs)
    ## hmm model
    essential_calls, _ = hmm_essential(disc_seq, regions)

# 2. Con-ARTIST part
## read data under experiment condition, the same as above
input_experiment = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\wig2\84_window_size=100.wig2'
experiment, *_ = SAMreader_Tn5(input_experiment, chr_search_term, genome_length,
                               chr_start, chr_end, chr_name, window_size=window_size)
if norm_flag:
    window_len = 100000
    normed = window_average(Tn5_NT, experiment, window_len, genome_length)
    experiment = list(normed)
## let all_Tn5_sum be the control
control = all_Tn5_sum
## Simulate bottleneck, set saturation_ratio=1 if no bottleneck is assumed.
num_boots = 100
saturation_ratio = (len(experiment) - experiment.count(0)) / (len(control) - control.count(0))
# saturation_ratio = 1
control_norm = simulate_equal_saturation(saturation_ratio, sum(experiment), control, num_boots)
## whether to use HMM
HMM_flag = True
p_value_threshold = 0.05
if HMM_flag:
    fc_thresholds = [-3, -1, 1, 3]
    min_sl, max_sl, min_enr, max_enr, prop_ignore = 0.5, 0.9, 0.5, 0.9, 0.1
    iter_flag = False  # whether iterate the HMMs, False in the original version.
    denominators = calc_cutoffs(experiment)
    state_confidence = hmm_train_mwu(experiment, control_norm, essential_calls, unique_indices, p_value_threshold,
                                     denominators, fc_thresholds, iter_flag=iter_flag)
    output_cond, call_stats = call_conditional(state_confidence, unique_indices, min_sl, max_sl, min_enr, max_enr,
                                               prop_ignore)
else:
    sig_proportion = 0.5
    mwu_boots = run_all_MWU_boots(control_norm, experiment, unique_indices)
    output_cond, call_stats = MWU_summary(mwu_boots, p_value_threshold, sig_proportion, unique_indices, control_norm,
                                          experiment)

# 3. output to tsv
file_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\ARTIST\Con-ARTIST\output_conditional_without_HMM_2.tsv'
kwargs = {'fi_path': [input_file, input_experiment],
          'samreader': [window_size],
          'norm': [norm_flag, window_len] if norm_flag else [norm_flag],
          'hmm_flag': HMM_flag,
          'p_value': p_value_threshold}
if HMM_flag:
    kwargs.update({
        'num_boots': num_boots,
        'disc_cond': [denominators, fc_thresholds],
        'iter_flag': iter_flag,
        'probs': [min_sl, max_sl, min_enr, max_enr, prop_ignore]})
else:
    kwargs.update({'sig_prop': sig_proportion})
output_to_tsv(file_path, unique_names, output_cond, call_stats, method='Con-ARTIST', **kwargs)
