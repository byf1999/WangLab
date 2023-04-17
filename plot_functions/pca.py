from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

input_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM'

df = pd.read_excel(os.path.join(input_path, r'PCA_input.xlsx')).values
df_norm = (df - np.mean(df, axis=0, keepdims=True)) / np.std(df, axis=0, keepdims=True)
# df_norm = (df - np.min(df, axis=1, keepdims=True)) / (np.max(df, axis=1, keepdims=True) - np.min(df, axis=1, keepdims=True))
# df_norm[np.isnan(df_norm)] = 0
pca = PCA(n_components=2)
o = pca.fit_transform(np.transpose(df_norm))
ratio = pca.explained_variance_ratio_

ax = plt.subplot()
ax.scatter(o[:3, 0], o[:3, 1], c='red')
ax.scatter(o[3:, 0], o[3:, 1], c='cyan')
ax.set_ylim(ax.get_xlim())
# ax.set_xticks([])
# ax.set_yticks([])
ax.set_xlabel(f'PC1: {ratio[0] * 100:.0f}% variance', size=12)
ax.set_ylabel(f'PC2: {ratio[1] * 100:.0f}% variance', size=12)
plt.show()

ax.get_figure().savefig(os.path.join(input_path, 'pics', 'pca.svg'))
