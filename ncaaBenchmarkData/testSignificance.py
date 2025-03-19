#!/dors/meilerlab/data/belle6/miniforge3/envs/glypred/bin/python

import matplotlib
matplotlib.rcParams["figure.dpi"] = 200
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import wilcoxon
from sys import argv

f=open(argv[1])
col1 = int(argv[2])
col2 = int(argv[3])
select = [int(n) for n in argv[2:]]

x = []
y = []
cols = []
labels = None
for line in f:
    if "ERROR" not in line:
        parts = line.split()
        if len(parts) <= 1:
            continue
        
        cols.append([float(i.split(":")[1]) if i.split(":")[1] != "" and i.split(":")[1] != "BAD" else -1. for i in parts[1:]])
        if labels == None:
            labels = [i.split(":")[0] for i in parts[1:]]
        poptime = False
        for colidx in select:
            if cols[-1][colidx] == -1.:
                cols.pop()
                break
        if poptime:
            continue
        try:
            x.append(float(cols[-1][col1]))
            y.append(float(cols[-1][col2]))
        except:
            print("Problem with %s"%parts[0])
print(wilcoxon(x,y))
print(sum(x)/len(x),sum(y)/len(y))
print(sum(x)/len(x)-sum(y)/len(y))
cols = np.array(cols)[:,select]
print(cols.shape)
q1, med, q3 = np.percentile(cols, [25, 50, 75], axis=1)

labels = [labels[i] for i in select]
fig = plt.figure()
ax = fig.add_subplot(111)
parts = ax.violinplot(cols,showmedians=True, showextrema=True, points=500)
for c in parts.keys():
    if c == "bodies":
        for pc in parts[c]:
            pc.set_alpha(1)
            pc.set_facecolor("orange")
            pc.set_edgecolor("black")
    else:
        parts[c].set_alpha(1)
        parts[c].set_color("black")
        if c == "cbars":
            parts[c].set_linestyle("dashed")
        if c == "cmins":
            parts[c].set_alpha(0)
        if c == "cmaxes":
            parts[c].set_alpha(0)

mode = argv[1].split(".")[0][-3:]
if mode == "prm" or mode == "rtm":
    ax.set_ylabel("Rotamer Recovery")
else:
    ax.set_ylabel("Sequence Recovery")
ax.set_ylim(bottom=0.0)
ax.set_xticks([i+1 for i in range(len(labels))])
ax.set_xticklabels(labels)
plt.show()

'''
l = []
for i in select:
    for j in range(cols.shape[0]):
        l.append(labels[i])

df = pd.DataFrame(np.array([cols.flatten(),l],dtype=object).T,columns=["Recovery","Method"])

df["Recovery"] = df["Recovery"].astype(float)
df["Method"] = df["Method"].astype(str)
print(df)
print(df.dtypes)
sns.violinplot(df=df,y="Recovery",x="Method")
'''
