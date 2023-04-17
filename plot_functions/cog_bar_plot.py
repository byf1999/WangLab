import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

fi_path = r'F:\LearningFiles\Master\8.TIS\A. ve 84筛选\SAM'
df = pd.read_table(os.path.join(fi_path, 'cog_input.txt'), sep='\t')
xticks = [annot[1:2] for annot in df['eggNOG_class_annotation']]

plot, (ax1, ax2) = plt.subplots(2,1, figsize=(3,7.5),sharex='all')
ax1.bar(xticks,df['genome_count'], color=sns.color_palette("hls", 22))
ax2.bar(xticks,df['E_count'], color=sns.color_palette("hls", 22))
plt.subplots_adjust(hspace=.05)
#ax1.set_xtickss
#ax1.spines['bottom'].set_visible(False)
#ax1.set_xticks(xticks,size=22)
ax2.set_xticks(xticks,size=22)
ax1.set_yticks(ax1.get_yticks(),size=22)
ax2.set_yticks(ax2.get_yticks(),size=22)
plt.show()
plot.savefig(os.path.join(fi_path, 'cog.svg'))