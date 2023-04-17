import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

fi_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM\FPKM.xlsx'
fn = pd.read_excel if os.path.splitext(fi_path)[1] in ['.xls', '.xlsx'] else pd.read_table

def plot(sheet_name=0):
    df = fn(fi_path, header=0, sheet_name=sheet_name)

    fc_threshold = 2
    ax = plt.subplot()
    ax.scatter(df.iloc[:, 2], df.iloc[:, 3], c=[190/256, 190/256, 190/256], s=5)
    ax.axline((np.log10(fc_threshold), 0), slope=1, c=[148/256,203/256,141/256], linestyle='--')
    ax.axline((0, np.log10(fc_threshold)), slope=1, c=[148/256,203/256,141/256], linestyle='--')
    ax.axline((0, 0), slope=1, c=[148/256,203/256,141/256], linestyle='--')
    ax.set_xlabel('log10 FPKM in control', size=12)
    ax.set_ylabel('log10 FPKM in experiment', size=12)
    ax.set_xlim([-1, 6])
    ax.set_ylim([-1, 6])
    ax.set_title(sheet_name)
    plt.show()

sheet_name = '26_M_End'
plot(sheet_name)
plot()