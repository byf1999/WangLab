import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

file_path = ''
df_fc = pd.read_excel(file_path)
df_norm = (df_fc - np.mean(df_fc, axis=0)) / np.std(df_fc, axis=0)

link_method, metric = 'average', 'mahalanobis'
threshold = 2.3
fig_size, dpi = (12.8, 7.2), 144
Z = linkage(df_norm, method=link_method, metric=metric, optimal_ordering=True)
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
'''
y_range = [-y_min / 2.1, -y_min / 21]
for i, flag in enumerate(df['t3']):
    if flag == True:
        ind = r['ivl'].index(str(i))
        x_range = [ind * 10 + 5] * 2
        ax.plot(x_range, y_range, 'r')
y_range = [-y_min, -y_min / 1.9]
for i, flag in enumerate(df['t6']):
    if flag == True:
        ind = r['ivl'].index(str(i))
        x_range = [ind * 10 + 5] * 2
        ax.plot(x_range, y_range, 'c')
'''