import argparse
params_dict = {
    ('input_path', 'i'): ('input file path', None),
    ('num_rep', 'nr'): ('number of replecation', None),
    ('num_point', 'np'): ('number of data point, either nr or np should be defined', None),
    ('output_path', 'o'): ('output file path', None),
    ('p_value', 'p'): ('threshold that assume a model is better than another, default 0.05', 0.05),
    ('confint', 'cf'): (
        'level of confidence interval, 0.95 means 95 percent of confidence interval, default 0.95', 0.95),
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
#parser = argparse.ArgumentParser(description='')
parser.add_argument('function', help='function used', default=None)
for (n1, n2), (h, default_value) in params_dict.items():
    parser.add_argument(f'--{n1}', f'-{n2}', help=f'{h}', default=default_value)
args = parser.parse_args()
