import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt


def set_params():
    params_dict = {
        ('input_path', 'i'): ('input file path', None),
        ('num_rep', 'nr'): ('number of replecation', None),
        ('num_point', 'np'): ('number of data point, either nr or np should be defined', None),
        ('output_path', 'o'): ('output file path', None),
        ('p_value', 'p'): ('threshold that assume a model is better than another, default 0.05', 0.05),
        ('confint', 'confint'): (
            'level of confidence interval, 0.95 means 95% confidence interval, default 0.95', 0.95),
        ('max_degree', 'max_degree'): ('max degree of linear model to be fitted, default 3', 3),
        ('link_method', 'link_method'): ('the linkage method used to compute the distance, default average', 'average'),
        ('gene_tag', 'tag'): ('gene tag that other genes are similar with', None),
        ('fc_threshold', 'fc_value'): ('threshold of l2fc, default -3', -3),
        ('dist_threshold', 'd_value'): ('threshold of distance between two genes, default 2', 2),
        ('cluster_threshold', 'cl_value'): ('threshold of clustering, default 3', 3),
        ('sample_norm_flag', 'sample_norm_flag'): (
            'flag whether the input data need normalization, if your input has been normalized, set this arg to False, default True',
            True),
        ('input_mode', 'mode'): (
            'the mode of input data, 0 means only 1 input column and 1 means input and output columns are one-to-one correspondent, default 0',
            0),
    }
    parser = argparse.ArgumentParser()
    # parser = argparse.ArgumentParser(description='')
    parser.add_argument('function', help='function used', default=None)
    for (n1, n2), (h, default_value) in params_dict.items():
        parser.add_argument(f'--{n1}', f'-{n2}', help=f'{h}', default=default_value)
    return parser.parse_args()


def FC_calulator(input, num_rep, filter=True, sample_norm_flag=True, mode=0):
    '''calculate average log2 fold change from the input file'''
    assert num_rep, 'nr should be defined!'

    df = pd.read_excel(input, header=0)
    data = df.values
    start = 2 if mode == 0 else 1 + num_rep
    col_output = data[:, start:].shape[1]
    times = np.unique([int(des.split('_')[0][1:]) for des in df.columns[start:]])
    times = np.array(times, dtype=np.int32).reshape((-1))  # time point
    if filter == True:
        # filter genes with insufficient data
        df['frac_data'] = np.sum(data[:, 2:] > 0, axis=1) / col_output
        data = df.loc[(df['frac_data'] > -1) & (df.iloc[:, 1] > -1)].values

    tag = data[:, 0]
    input = (data[:, 1:start]).reshape((-1, start - 1))
    output = data[:, start:(col_output + start)]
    # sample normalization
    if sample_norm_flag:
        input = input / np.sum(input, axis=0, keepdims=True) * 1e8
        output = output / np.sum(output, axis=0, keepdims=True) * 1e8

    # input = np.repeat(input, col_output/num_rep).reshape(output.shape)
    input = np.repeat(input, col_output / (start - 1)).reshape(output.shape)
    l2fc = np.log2(((output + 1) / (input + 1)).astype('float'))
    if num_rep > 1:
        mean_fc = []
        sd_fc = []
        for i in range(0, col_output, num_rep):
            mean_fc.append(np.mean(l2fc[:, i:i + num_rep], axis=1))
            sd_fc.append(np.std(l2fc[:, i:i + num_rep], axis=1, ddof=1))
        mean_fc = np.array(mean_fc)
        sd_fc = np.array(sd_fc)
        sd_fc[sd_fc == 0] = np.max(sd_fc)  # 有待商榷
        w = 1 / (sd_fc ** 2)
    elif num_rep == 1:
        mean_fc = np.transpose(l2fc)
        w = np.ones_like(mean_fc)
    else:
        raise ValueError('-nr should not be less than 1')
    return (times, tag, mean_fc, w)


def lm(x, y, w, degree):
    '''polynomial model fitting (y ~ x ^ degree) using weighted least-square'''
    model = sm.WLS(y, np.vander(x, degree + 1, increasing=True), weights=w)
    res = model.fit()
    return res


def model_fit(mean_fc, times, w, max_degree, p_value, confint):
    '''
    fit the data using polynomial model with maximum model degree of max_degree
    :param mean_fc: Y
    :param times: X
    :param w: weight
    :param max_degree: maximum model degree
    :param p_value: p_value threshold to define a model is better than another
    :param confint: percentage of confident interval
    :return: best model's parameters and the upper confident interval of slope(b, where y ~ a + bx + cx^2 + ...)
    '''
    num_genes = mean_fc.shape[1]
    max_degree = min(max_degree, mean_fc.shape[0] - 1)
    best_res = [lm(times, mean_fc[:, i], w[:, i], degree=0) for i in range(num_genes)]
    best_degrees = np.zeros(num_genes)
    slopeUpperCI = []
    best_params = []
    for i in range(num_genes):
        for degree in range(1, max_degree + 1):
            res = lm(times, mean_fc[:, i], w[:, i], degree=degree)
            if res.compare_f_test(best_res[i])[1] < p_value:  ###
                best_res[i] = res
                best_degrees[i] = degree
        if best_degrees[i] == 0:
            slopeCI = [0, 0]
        else:
            slopeCI = best_res[i].conf_int(alpha=1 - confint)[1,]
        slopeUpperCI.append(slopeCI[1])
        best_params.append(best_res[i].params)
        for _ in range(max_degree - int(best_degrees[i])):
            best_params[i] = np.append(best_params[i], 0)
    return (best_params, slopeUpperCI)


def plot(Z, df, threshold=2.3, fig_size=(12.8, 7.2), dpi=144):
    '''plot dendrogram of hierarchical clustering results'''
    from scipy.cluster.hierarchy import dendrogram
    plt.figure(figsize=fig_size, dpi=dpi)
    r = dendrogram(Z, color_threshold=threshold)

    x_min = -Z.shape[0] / 10
    x_max = Z.shape[0] * 10 + 10
    y_max = np.ceil(np.max(r['dcoord']))
    y_min = y_max / 20
    plt.xlim(x_min, x_max);
    plt.ylim(-y_min, y_max);
    plt.xlabel('');
    plt.ylabel('Phenotypic distance');
    plt.xticks([], []);
    # plt.yticks([])
    ax = plt.gca()
    # for s in ['top', 'bottom', 'right', 'left']:
    for s in ['top', 'bottom', 'right']:
        ax.spines[s].set_visible(False)

    # plot t3 and t6 genes
    y_range = [-y_min / 1.9, -y_min / 21]
    for i, flag in enumerate(df['t3']):
        if flag == True:
            ind = r['ivl'].index(str(i))
            x_range = [ind * 10 + 5] * 2
            ax.plot(x_range, y_range, 'r')
    y_range = [-y_min, -y_min / 2.1]
    for i, flag in enumerate(df['t6']):
        if flag == True:
            ind = r['ivl'].index(str(i))
            x_range = [ind * 10 + 5] * 2
            ax.plot(x_range, y_range, 'c')
