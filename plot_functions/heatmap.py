import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, palettable
import seaborn as sns
import seaborn.matrix

input_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM'

df = pd.read_excel(os.path.join(input_path, r'PCA_input.xlsx')).values
df_norm = (df - np.mean(df, axis=0, keepdims=True)) / np.std(df, axis=0, keepdims=True)
r = np.corrcoef(np.transpose(df_norm))
r = np.corrcoef(np.transpose(df))
df = pd.DataFrame(r)
groups = ['c1', 'c2', 'c3', 'e1', 'e2', 'e3']
df.columns = groups
df.index = groups
# cmaps = ['Blues', 'Blues_r', 'OrRd', 'OrRd_r', 'PuBu', 'PuBu_r', ]

plt.clf()
plt.figure(dpi=120)
ax = sns.heatmap(data=df, cmap='OrRd', annot=True, fmt='.2f', cbar=True,
                 cbar_kws={  # 'label': 'ColorbarName',  # color bar的名称
                     # 'orientation': 'horizontal',  # color bar的方向设置，默认为'vertical'，可水平显示'horizontal'
                     "ticks": np.arange(0.85, 1.0, 0.05),  # color bar中刻度值范围和间隔
                     "format": "%.2f",  # 格式化输出color bar中刻度值
                     "pad": 0.1,  # color bar与热图之间距离，距离变大热图会被压缩
                 },
                 )
# ax.imshow(df, cmap='OrRd')
plt.show()
plot = ax.get_figure()
plot.savefig('test.svg')

plt.figure(dpi=200)
# sns.clustermap(pd.DataFrame(np.transpose(df_norm)))
# sns.clustermap(df_norm)
ax = sns.clustermap(df, cmap='OrRd', figsize=(10, 10),
                    # cbar_pos=(0.85, 0.35, 0.03, 0.45),  # (left, bottom, width, height))
                    cbar_pos=(0.06, 0.7, 0.03, 0.25),  # (left, bottom, width, height))
                    tree_kws={'linestyles': 'solid',  # 线型
                              'colors': 'black',  # 线色
                              'linewidths': 1},  # 线宽
                    annot=True, fmt='.2f', cbar=True,
                    cbar_kws={  # 'label': 'ColorbarName',  # color bar的名称
                        # 'orientation': 'horizontal',  # color bar的方向设置，默认为'vertical'，可水平显示'horizontal'
                        "ticks": np.arange(0.85, 1.0, 0.05),  # color bar中刻度值范围和间隔
                        "format": "%.2f",  # 格式化输出color bar中刻度值
                        "pad": 0.3,  # color bar与热图之间距离，距离变大热图会被压缩
                    },
                    annot_kws={'fontsize': 20},
                    xticklabels=groups)
# plot1 = ax.figure
# plot1.xticks(fontsize=20)
plt.show()
plot.savefig('test1.svg')
ax.savefig('test1.svg')
