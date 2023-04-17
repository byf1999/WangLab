'''
Dependency: numpy, pandas, scikit-learn
'''
import re
import time

from sklearn.neighbors import KernelDensity
import numpy as np
import pandas as pd
import collections
import random
from scipy.stats import ranksums


# from scipy.stats import norm, gaussian_kde
# from scipy.special import ndtr


def genome_parser(gtf_file, chr_search_term, chrom_end, cds_search_term=r'ID=(.+?);'):
    df = pd.read_table(gtf_file, sep='\t', header=None)
    df = df.loc[(df.iloc[:, 0] == chr_search_term) & (df.iloc[:, 2] == 'CDS')]
    df = df.sort_values(by=3)  # sort by start pos

    df[8] = [re.search(cds_search_term, des).group(1) for des in df.iloc[:, 8]]
    IGname, IGstart, IGend = ['IG_1'], [1], [0]
    for i in range(len(df) - 1):
        if df.iloc[i + 1, 3] > df.iloc[i, 4] + 1:
            IGname.append('IG' + df.iloc[i + 1, 8])
            IGstart.append(df.iloc[i, 4] + 1)
            IGend.append(df.iloc[i + 1, 3] - 1)
    if df.iloc[-1, 4] < chrom_end - 1:
        IGname.append('IG_chrm_end')
        IGstart.append(df.iloc[-1, 4] + 1)
        IGend.append(chrom_end)
    df1 = pd.DataFrame({3: IGstart, 4: IGend, 8: IGname})
    df = pd.concat([df.iloc[:, [3, 4, 8]], df1])

    df = df.sort_values(by=[3, 4])
    return (list(df.iloc[:, 0]), list(df.iloc[:, 1]), list(df.iloc[:, 2]))


def SAMreader_Tn5(sam_file, chr_identifier, genome_length, chr_start, chr_end, chr_name, window_size=100):
    # generate Tn5_annotations and Tn5_NT
    Tn5_NT = list(range(round(window_size / 2), genome_length, window_size))
    if Tn5_NT[-1] + round(window_size / 2)  < genome_length:
        Tn5_NT.append(genome_length)
    Tn5_annotations = []
    s = 0
    for nt in Tn5_NT:
        for i in range(s, len(chr_name)):
            if chr_start[i] <= nt and chr_end[i] >= nt:
                Tn5_annotations.append(chr_name[i])
                s = i
                break

    # calculate all_Tn5_sum
    if sam_file.endswith('.wig2'):  # if input a .wig2 file
        all_Tn5_sum = list(pd.read_table(sam_file, header=None).values.squeeze())
        return (all_Tn5_sum, Tn5_NT, Tn5_annotations)
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

    all_Tn5_sum = [sum(all_pos0[i - 1:i - 1 + window_size]) for i in range(1, genome_length + 1, window_size)]
    return (all_Tn5_sum, Tn5_NT, Tn5_annotations)


def get_unique_names(Tn5_annotations):
    unique_annotations = list(set(Tn5_annotations))
    unique_annotations.sort(key=Tn5_annotations.index)
    counts = [Tn5_annotations.count(annotation) for annotation in unique_annotations]
    unique_indices = np.cumsum(counts)
    tmp = 1 + unique_indices
    unique_indices = [[1] + list(tmp[:-1]), list(unique_indices)]
    return (unique_annotations, unique_indices)


def window_average(Tn5_NT, all_Tn5_sum, window_len, genome_length):
    chr_avg = sum(all_Tn5_sum) / (len(all_Tn5_sum) - all_Tn5_sum.count(0))
    window_interval = list(range(0, genome_length, window_len))
    if window_interval[-1] < genome_length: window_interval += [genome_length]
    window_indices = []
    window_avg = []
    ind = 0
    s, count = 0, 0
    for i, nt in enumerate(Tn5_NT):
        if nt > window_interval[ind + 1]:
            ind += 1
            window_avg.append(s / count)
            s, count = 0, 0
        window_indices.append(ind)
        s += all_Tn5_sum[i]
        count += 1 if all_Tn5_sum[i] else 0
    window_avg.append(s / count)
    window_scale = chr_avg / np.array(window_avg)
    scale = [window_scale[ind] for ind in window_indices]
    return np.array(np.round(np.multiply(scale, all_Tn5_sum)), dtype='int')


def sliding_window(all_Tn5_sum, p_value=0.03, sim_times=1000, sw_size=7):
    # bootstrap
    total_sum = [sum(random.sample(all_Tn5_sum, sw_size)) for _ in range(sim_times)]

    # assign p_value of each window
    essential_regions = np.zeros((len(all_Tn5_sum), 1))
    essential_pvalues = np.zeros((len(all_Tn5_sum), 1))
    # kde = gaussian_kde(total_sum)
    kde = KernelDensity(bandwidth=0.5).fit(np.array(total_sum).reshape((-1, 1)))
    cdf = {}
    i = 0
    h = 0.01
    while True:
        cdf[i] = sum(np.exp(kde.score_samples(np.arange(-2, i, h).reshape((-1, 1))))) * h
        if cdf[i] > p_value:
            break
        i += 1
    m = cdf[i]
    for i in range(0, len(all_Tn5_sum) - sw_size + 1):
        s = int(np.sum(all_Tn5_sum[i:i + sw_size]))
        if cdf.get(s, m) < p_value:
            essential_regions[i:i + sw_size] = 1
        essential_pvalues[i:i + sw_size] = cdf.get(s, m)

    l = 1
    for i in range(1, sw_size):
        s = int(np.sum(all_Tn5_sum[-i:]) * sw_size / l)
        l += 1
        if cdf.get(s, m) < p_value:
            essential_regions[-i:] = 1
        essential_pvalues[-i:] = cdf.get(s, m)
    return essential_regions, essential_pvalues


def calc_cutoffs(x, whisker=1.5):
    c25, c75 = np.percentile(x, [25, 75])
    up_outlier = c75 + whisker * (c75 - c25)
    if c25 == 0:
        c25 = 1
    if c75 == c25:
        c75 += 1
    if up_outlier == c75:
        up_outlier += 1
    return [c25, c75, up_outlier]


def discretize(all_Tn5_sum, zero, c25, c75, c_out):
    disc_data = np.zeros_like(all_Tn5_sum)
    for i, num in enumerate(all_Tn5_sum):
        if num <= zero:
            continue
        elif num <= c25:
            disc_data[i] = 1
        elif num <= c75:
            disc_data[i] = 2
        elif num <= c_out:
            disc_data[i] = 3
        else:
            disc_data[i] = 4
    return disc_data


def hmm_estimate(disc_seq, state, num_state=0, num_symbol=0):
    num_state = num_state if num_state else int(np.max(state)) + 1
    num_symbol = num_symbol if num_symbol else int(np.max(disc_seq)) + 1

    tr = np.zeros((num_state, num_state))
    E = np.zeros((num_state, num_symbol))

    for i in range(len(disc_seq) - 1):
        tr[int(state[i]), int(state[i + 1])] += 1
        E[int(state[i]), int(disc_seq[i])] += 1
    E[int(state[i + 1]), int(disc_seq[i + 1])] += 1

    tr = tr / np.sum(tr, axis=1, keepdims=True)
    E = E / np.sum(E, axis=1, keepdims=True)
    tr[np.isnan(tr)] = 0
    E[np.isnan(E)] = 0
    return (tr, E)

    # state = [1 if i==2 else 0 for i in [1,1,1,1,2,2,2,1,1,2]]
    # disc_seq = list(np.array([1,1,1,2,1,2,2,3,4,5])-1)


def hmm_essential(disc_seq, state):
    tr, E = hmm_estimate(disc_seq, state)
    state, count = hmm_converge(disc_seq, tr, E)
    return state, count


def hmm_converge(disc_seq, tr, E):
    a0 = np.zeros_like(tr)
    count = 0
    state = None
    while np.sum(a0 - tr) != 0:
        a0 = tr
        state = hmm_viterbi(disc_seq, tr, E)
        tr, E = hmm_estimate(disc_seq, state, *E.shape)
        count += 1
    return state, count


def hmm_viterbi(disc_seq, tr, E):
    num_state = tr.shape[0]
    l = len(disc_seq)
    logTR = np.log(tr)
    logE = np.log(E)
    pTR = np.zeros((num_state, l))
    v = [-np.inf] * num_state
    v[0] = 0
    vOld = v.copy()
    for count in range(l):
        for state in range(num_state):
            bestVal = -np.inf
            bestPTR = 0
            for inner in range(num_state):
                val = vOld[inner] + logTR[inner, state]
                if val > bestVal:
                    bestVal = val
                    bestPTR = inner
            pTR[state, count] = bestPTR
            v[state] = logE[state, disc_seq[count]] + bestVal
        vOld = v.copy()
    final_state = np.argmax(v)

    current_state = np.zeros(l, dtype='int')
    current_state[-1] = final_state
    for i in range(2, l):
        current_state[-i] = pTR[current_state[-i + 1], -i + 1]
    return current_state


def output_in_cds(essential_calls, unique_indices):
    unique_calls = np.zeros_like(unique_indices[0])
    for i, (start, end) in enumerate(zip(unique_indices[0], unique_indices[1])):
        calls = essential_calls[start - 1:end]
        l = len(calls)
        mid_calls = calls[round(l / 5):round(l * 0.8)]
        if any(mid_calls): unique_calls[i] = 2
        if all(mid_calls): unique_calls[i] = 1
    unique_calls = list(unique_calls)
    out_stats = [unique_calls.count(i) for i in range(3)]
    return unique_calls, out_stats


def write_el_artist(file_path, unique_names, unique_calls, out_stats, fi_path, samreader, norm, sw, dics, sep):
    with open(file_path, 'w') as fo:
        fo.write('# Using EL-ARTIST pipeline.\n')
        fo.write(f'# Input file: {fi_path[0]}\n')
        fo.write(f'# Parameters: window size of SAMreader: {samreader[0]} bp, with saturation of {samreader[1]}.\n')
        fo.write(f'# Normalizaion: {norm[0]}' + (f', with window_len of {norm[1]} bp' if norm[0] else '') + '\n')
        fo.write(f'# Sliding window: simulation times {sw[1]}, window size {sw[0]}, and p value threshold {sw[2]}\n')
        fo.write(f'# Discritize: with cutoffs {dics[0]}\n')
        fo.write(
            f'# essential regions: {out_stats[1]}, domain-essential: {out_stats[2]}, non-essential: {out_stats[0]}\n')
        fo.write(f'gene_id{sep}call\n')
        for name, call in zip(unique_names, unique_calls):
            fo.write(f'{name}{sep}{call}\n')


def write_con_artist(file_path, unique_names, unique_calls, out_stats, fi_path, samreader, norm, sw, dics, hmm_flag,
                     num_boots, disc_cond, p_value, iter_flag, probs, sig_prop, sep):
    with open(file_path, 'w') as fo:
        fo.write('# Using Con-ARTIST pipeline.\n')
        fo.write(f'# Control file: {fi_path[0]};\n')
        fo.write(f'# Experiment file: {fi_path[1]};\n')
        fo.write(f'# Parameters of Con-ARTIST:\n')
        fo.write(f'## Window size of SAMreader: {samreader[0]} bp.\n')
        fo.write(f'## Normalizaion: {norm[0]}' + (f', with window_len of {norm[1]} bp' if norm[0] else '') + '\n')
        fo.write(f'## HMM used: {hmm_flag}\n')
        fo.write(f'## P-value threshold: {p_value}\n')
        if hmm_flag:
            headers = ['gene_id', 'calls']
            type_names = ['domain conditionally essential',
                          'conditionally essential',
                          'domain conditionally enriched',
                          'conditionally enriched',
                          'not different']
            fo.write(f'## times of simulating bottleneck: {num_boots}\n')
            fo.write(f'## Discretize: with cutoffs {disc_cond[0]} and fc_cutoffs {disc_cond[1]}\n')
            fo.write(f'## Iteration during HMM: {iter_flag}\n')
            fo.write(f'## Probability to decide: {probs}\n')
        else:
            headers = ['gene_id', 'percent of calls', 'calls',
                       'mean p_value', 'std p_value',
                       'mean fold change', 'std fold change']
            type_names = ['not different', 'significant different']
            fo.write(f'## Significant proportion is: {sig_prop}\n')
        for i, t in enumerate(type_names):
            fo.write(f'# {t}: {out_stats[i]}\n')
        fo.write(sep.join(headers) + '\n')
        for name, call in zip(unique_names, np.reshape(unique_calls, (len(unique_calls), -1))):
            line = f'{name}'
            for c in call:
                line += f'{sep}{c}'
            fo.write(f'{line}\n')


def output_to_tsv(file_path, unique_names, unique_calls, out_stats, method, sep='\t', **kwargs):
    samreader = kwargs.get('samreader')
    norm = kwargs.get('norm')
    sw = kwargs.get('sw')
    dics = kwargs.get('discritize')
    fi_path = kwargs.get('fi_path')
    if method.upper() == 'EL-ARTIST':
        write_el_artist(file_path, unique_names, unique_calls, out_stats, fi_path, samreader, norm, sw, dics, sep)
    elif method.upper() == 'CON-ARTIST':
        hmm_flag = kwargs.get('hmm_flag')
        num_boots = kwargs.get('num_boots')
        disc_cond = kwargs.get('disc_cond')
        p_value = kwargs.get('p_value')
        iter_flag = kwargs.get('iter_flag')
        probs = kwargs.get('probs')
        sig_prop = kwargs.get('sig_prop')
        write_con_artist(file_path, unique_names, unique_calls, out_stats, fi_path, samreader, norm, sw, dics, hmm_flag,
                         num_boots, disc_cond, p_value, iter_flag, probs, sig_prop, sep)


def simulate_equal_saturation(saturation_ratio, num_reads_in_exp, control, num_boots):
    assert saturation_ratio <= 1, 'the ratio of saturation between experiment and control group is less then 1, please check.'
    input_proportion_norm = control / sum(control) * saturation_ratio
    input_proportion_norm = np.concatenate([input_proportion_norm, np.array([1 - saturation_ratio])], axis=0)
    input_samples = np.random.multinomial(num_reads_in_exp, input_proportion_norm, size=num_boots)
    weights = num_reads_in_exp / (num_reads_in_exp - input_samples[:, -1])
    return np.multiply(np.transpose(input_samples), weights)[:-1, ]


def run_all_MWU_boots(control_norm, experiment, unique_indices):
    p_values = []
    for i in range(np.shape(control_norm)[1]):
        p_values.append(MWU_test(experiment, control_norm[:, i], unique_indices))
    return np.transpose(np.array(p_values))


def hmm_train_mwu(experiment, control_norm, essential_calls, unique_indices, mwu_threshold, denominators,
                  fc_thresholds, iter_flag=False):
    state_boots = []
    for i in range(np.shape(control_norm)[1]):
        control = control_norm[:, i]
        disc_seq = discretize_conditional(experiment, control, denominators, fc_thresholds)
        init_state = essential_calls
        for start, end in zip(unique_indices[0], unique_indices[1]):
            score_1 = experiment[start - 1:end]
            score_2 = control[start - 1:end]
            _, p = ranksums(score_1, score_2)
            if p < mwu_threshold:
                if sum(score_1) > sum(score_2):
                    init_state[start - 1:end] = 2
                elif sum(score_1) < sum(score_2):
                    init_state[start - 1:end] = 3
        tr, E = hmm_estimate(disc_seq, init_state, num_state=4, num_symbol=15)
        if iter_flag:
            state, _ = hmm_converge(disc_seq, tr, E)
        else:
            state = hmm_viterbi(disc_seq, tr, E)
        state_boots.append(state)
    state_boots = np.transpose(np.array(state_boots))
    state_confidence = []
    for i in range(state_boots.shape[0]):
        tmp = collections.Counter({i: 0 for i in range(4)})
        tmp.update(collections.Counter(state_boots[i, :]))
        state_confidence.append(list(tmp.values()))
    state_confidence = np.array(state_confidence)
    return state_confidence / np.sum(state_confidence, axis=1, keepdims=True)


def call_conditional(state_confidence, unique_indices, min_sl, max_sl, min_enr, max_enr, prop_ignore):
    output_cond = []
    for start, end in zip(unique_indices[0], unique_indices[1]):
        call_sl = state_confidence[start - 1:end, 3]
        call_enr = state_confidence[start - 1:end, 2]
        l = end - start + 1
        call_sl = call_sl[round(prop_ignore * l): round((1 - prop_ignore) * l)]
        if all(call_sl > min_sl):
            cond = 1
        elif any(call_sl > max_sl):
            cond = 0
        elif all(call_enr > min_enr):
            cond = 3
        elif any(call_enr > max_enr):
            cond = 2
        else:
            cond = 4
        output_cond.append(cond)
    call_stats = collections.Counter(output_cond)
    call_stats.update({i: 0 for i in range(5)})
    return output_cond, call_stats


def discretize_conditional(experiment, control, denominators, fc_thresholds):
    denominator1, denominator2, denominator3 = denominators
    q1, q2, q3, q4 = fc_thresholds
    disc_seq = np.zeros_like(control, dtype=int)
    for i in range(np.shape(control)[0]):
        exp = experiment[i]
        ctrl = control[i]
        if ctrl:  # control is not 0
            l2fc = np.log2(exp / ctrl)
            if l2fc < q1:
                disc_seq[i] = 6
            elif l2fc < q2:
                disc_seq[i] = 5
            elif l2fc < 0:
                disc_seq[i] = 4
            elif l2fc == 0:  # experiment is 0
                if ctrl < denominator1:
                    disc_seq[i] = 11
                elif ctrl < denominator2:
                    disc_seq[i] = 12
                elif ctrl < denominator3:
                    disc_seq[i] = 13
                else:
                    disc_seq[i] = 14
            elif l2fc < q3:
                disc_seq[i] = 1
            elif l2fc < q4:
                disc_seq[i] = 2
            else:
                disc_seq[i] = 3
        else:  # control is 0
            if exp > denominator3:
                disc_seq[i] = 10
            elif exp > denominator2:
                disc_seq[i] = 9
            elif exp > denominator1:
                disc_seq[i] = 8
            elif exp > 0:
                disc_seq[i] = 7
    return disc_seq


def MWU_test(experiment, control, unique_indices):
    p_values = []
    for start, end in zip(unique_indices[0], unique_indices[1]):
        score_1 = experiment[start - 1:end]
        score_2 = control[start - 1:end]
        _, p = ranksums(score_1, score_2)
        p_values.append(p)
    return p_values


def MWU_summary(mwu_boots, p_value_threshold, sig_prop, unique_indices, control_norm, experiment):
    summary = np.zeros((np.shape(mwu_boots)[0], 6))
    summary[:, 0] = np.mean(mwu_boots < p_value_threshold, axis=1)
    summary[:, 1] = np.array(summary[:, 0] > sig_prop, dtype=int)
    summary[:, 2] = np.mean(mwu_boots, axis=1)
    summary[:, 3] = np.std(mwu_boots, axis=1)

    ratios = np.zeros_like(mwu_boots)
    for i, (start, end) in enumerate(zip(unique_indices[0], unique_indices[1])):
        exp_sum = sum(experiment[start - 1:end])
        for j in range(np.shape(control_norm)[1]):
            control = control_norm[:, j]
            ctrl_sum = sum(control[start - 1:end])
            if exp_sum == 0:
                if ctrl_sum == 0:
                    ratio = 0
                else:
                    ratio = 1 / ctrl_sum
            elif ctrl_sum == 0:
                ratio = exp_sum
            else:
                ratio = exp_sum / ctrl_sum
            ratios[i, j] = ratio

    summary[:, 4] = np.mean(ratios, axis=1)
    summary[:, 5] = np.std(ratios, axis=1)
    c = collections.Counter({0: 0, 1: 0})
    c.update(list(summary[:, 1]))
    return summary, c


def load_essential_calls(file_path):
    essential_calls = []
    with open(file_path, 'r') as fi:
        while True:
            content = fi.readline()
            if not content:
                break
            if content.startswith('#'): continue
            call = content.split('\t')[-1]
            essential_calls.append(call)
    return np.array([int(call) for call in essential_calls[1:]])
