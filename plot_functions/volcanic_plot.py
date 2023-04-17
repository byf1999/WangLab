import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

fi_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM'
file_name = r'volcanic_input.xlsx'
fn = pd.read_excel if os.path.splitext(fi_path)[1] in ['.xls', '.xlsx'] else pd.read_table

df = fn(os.path.join(fi_path, file_name), header=0)
df['-lgp'] = -np.log10(df.iloc[:,1])
p_threshold = 0.05
l2fc_threshold = 1

enr = (df.iloc[:,0] >= l2fc_threshold) & (df.iloc[:,1] < p_threshold)
dec = (df.iloc[:,0] <= -l2fc_threshold) & (df.iloc[:,1] < p_threshold)

ax = plt.subplot()
ax.scatter(df.iloc[:,0], df.iloc[:,2], c='grey',s=5)
ax.scatter(df[enr].iloc[:, 0], df[enr].iloc[:, 2], c='red', s=5)
ax.scatter(df[dec].iloc[:, 0], df[dec].iloc[:, 2], c='cyan',s=5)
a, b = ax.get_xlim()
x_max = max(abs(a), b)
ax.set_xlim([-x_max, x_max])
ax.set_xlabel('L2FC');
ax.set_ylabel('-log10 p_value');
ax.set_title('Volcanic plot')
ax.axvline(-1, linestyle='--');
ax.axvline(1, linestyle='--');
ax.axhline(-np.log10(0.05), linestyle='--');
plt.show()

ax.get_figure().savefig(os.path.join(fi_path, 'pics', 'volcanic_plot.svg'))