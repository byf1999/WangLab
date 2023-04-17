import matplotlib.pyplot as plt
import pandas as pd

fi_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\wig2\WT_window_size=100.wig2'
df = pd.read_table(fi_path,header=None)
ax = plt.subplot()
ax.scatter(df.index*100, df[0])
ax.set_ylim(0, 1e4)
ax.set_xlabel('Gene locus (bp)',size=12)
ax.set_ylabel('# of insertions',size=12)
plt.show()
ax.get_figure().savefig(r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM\pics\scatter.svg')
