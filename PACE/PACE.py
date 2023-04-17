'''
PACE python 3.x version

usage: python PACE.py function_name -i input_path.xlsx -nr 3 [function params]

input_path.xlsx format: each row represent a gene and each column represent a sample,
first row is column names, fist column is gene name and the second is input_path reads,
the following sample columns should be named like D1_OrganName_ParaName.
example:
        tag     input   D1_liver_a  D1_liver_b  D2_liver_a   D2_liver_b ...
        gene1   20      10          12          8            7
        gene2   100     50          48          20           18
        ...

package dependency: numpy, pandas, statsmodels, scipy, matplotlib,

example:
    python PACE.py IVD-finder -i input.xlsx -nr 3 -o output
    python PACE.py SimGenFinder -i input.xlsx -nr 3 -tag ETAE_RS11320
    python PACE.py Cluster -i input.xlsx -nr 3 -cl_value 2

'''

from utils import *
import os
import time

args = set_params()
input_path = args.input_path
output_path = args.output_path
num_rep = int(args.num_rep)
p_value = float(args.p_value)
confint = float(args.confint)
max_degree = int(args.max_degree)
link_method = args.link_method
cluster_threshold = float(args.cluster_threshold)
sample_norm_flag = args.sample_norm_flag
input_mode = args.input_mode


def IVD_finder(df):
    '''Restrict attention to IVD genes'''
    # 可进一步分类，而不是一刀切
    cut_index = df['slopeUpperCI'] < 0
    df = df[cut_index]
    return df, cut_index


def similar_gene_finder(best_params, tag, mean_fc, gene_tag, metric='mahalanobis', ivd_cut=None):
    '''find genes similar to a specified gene based on the distance calculated by L2FC'''
    from scipy.spatial.distance import cdist
    fc_threshold = args.fc_threshold
    dist_threshold = args.dist_threshold
    distance = cdist(np.array(best_params), np.array(best_params), metric=metric)
    index = list(tag).index(gene_tag)
    d_gene = distance[index]
    df = pd.DataFrame({'tag': tag, 'distance': d_gene, 'L2FC': np.mean(mean_fc[2:], axis=0)})
    if ivd_cut:
        df1 = df[ivd_cut]
        df2 = df[(df['L2FC'] < fc_threshold) & (df['distance'] < dist_threshold)]
        plt.scatter(df1['L2FC'], df1['distance'], color='grey');
        plt.scatter(df2['L2FC'], df2['distance'], color='r');
        a, b = np.where(df == tag)
        plt.scatter(df.iloc[int(a),]['L2FC'], df.iloc[int(a),]['distance'], color='c');
        plt.xlabel('Outbreak Index (log2(FC))');
        plt.ylabel('Distance from aroC profile(SD)');
        plt.vlines(fc_threshold, plt.ylim()[0], plt.ylim()[1], linestyle='dashed');
        plt.hlines(dist_threshold, plt.xlim()[0], plt.xlim()[1], linestyle='dashed');
        plt.savefig(os.path.join(os.getcwd(), fr'output\similar_gene.svg'))
    cut_index = (df['distance'] < dist_threshold) & (df['L2FC'] < fc_threshold)
    df = df[cut_index]
    return df


def gene_cluster(df, link_method='average', metric='mahalanobis', threshold=2.9):
    '''hierarchical clustering based on model params'''
    from scipy.cluster.hierarchy import linkage, fcluster
    dta = df.iloc[:, :(max_degree + 1)]
    df_norm = (dta - np.mean(dta, axis=0)) / np.std(dta, axis=0)
    df['tag_num'] = [int(t[-5:]) for t in df['tag']]
    df['t3'] = (df['tag_num'] <= 4215) & (df['tag_num'] >= 4070)
    df['t6'] = (df['tag_num'] <= 11310) & (df['tag_num'] >= 11235)
    Z = linkage(df_norm, method=link_method, metric=metric, optimal_ordering=True)
    plot(Z, df, threshold)
    # print(flag)
    while False:
        flag = input('Enter a threshold to generate a new graph, or enter y to finish.')
        if flag == 'y':
            break
        threshold = float(flag)
        plot(Z, df, threshold)
    plt.savefig(os.path.join(os.getcwd(), fr'output\cluster.svg'))
    cr = fcluster(Z, threshold, criterion='distance')
    df['cr'] = cr
    df = df.drop('tag_num', axis=1)
    # print(max(df['cr']))
    return df


def main(max_degree):
    # input_path = r'F:\LearningFiles\Master\18.LearningMaterial\Python\PACE\input.xlsx'
    # num_rep, p_value, confint, max_degree = 3, .05, 0.95, 3
    # link_method, metric = 'average', 'mahalanobis'
    # sample_norm_flag = True
    print(f'Running {args.function}......')
    t0 = time.time()
    norm_flag = False
    assert input_path is not None, 'Please enter input file path!'
    times, tag, mean_fc, w = FC_calulator(input_path, num_rep, sample_norm_flag=sample_norm_flag,
                                          mode=input_mode)
    raw_df = pd.DataFrame(mean_fc.transpose())
    raw_df['tag'] = tag

    if norm_flag:
        def normalization(mean_fc):
            data = mean_fc - np.mean(mean_fc[:, np.max(abs(mean_fc), axis=0) < 2], axis=1, keepdims=True)
            for i, index in enumerate(range(0, data.shape[1], 30)):
                ax = plt.subplot(10, 11, i + 1)
                ax.plot(data[:, index:index + 30])
            return data

        mean_fc = normalization(mean_fc)
    # Fit with linear model
    ## model params and slopeUpperCI are needed
    best_params, slopeUpperCI = model_fit(mean_fc, times, w, max_degree=max_degree, p_value=p_value, confint=confint)

    df = pd.DataFrame(np.array(best_params), columns=['intercept', 'slope', 'param3', 'param4'][:(max_degree + 1)])
    df['tag'] = tag
    df['slopeUpperCI'] = slopeUpperCI
    df, cut = IVD_finder(df)

    file_path = os.path.join(os.getcwd(), r'output')
    if not os.path.exists(file_path):
        os.mkdir(file_path)

    if args.function == 'IVD-finder':
        pass
    elif args.function == 'SimGenFinder':
        assert args.gene_tag is not None, 'Please enter a gene tag'
        df = similar_gene_finder(best_params, tag, mean_fc, args.gene_tag, ivd_cut=cut)
    elif args.function == 'Cluster':
        df = gene_cluster(df, link_method=link_method, threshold=cluster_threshold)
    else:
        print('Function not defined, please check your command!')
        return

    if output_path is None:
        output_name = args.function
    else:
        output_name = output_path

    df.to_excel(os.path.join(os.getcwd(), fr'output\{output_name}.xlsx'))
    print(f'Task finished, cost {time.time() - t0} s')
    raw_df.to_excel(os.path.join(os.getcwd(), fr'output\raw_L2FC.xlsx'))


if __name__ == '__main__':
    main(max_degree)
