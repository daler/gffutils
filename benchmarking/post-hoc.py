import pandas
import seaborn as sns
import sqlite3
import json
from matplotlib import pyplot as plt
import numpy as np

df = pandas.read_sql_query('select * from timings', sqlite3.connect('timings.db'))

dd = []
for _, row in df.iterrows():
    d = row.to_dict()
    k = json.loads(row.kwargs)
    for key, val in k.items():
        d[key] = val
    for key, val in k['pragmas'].items():
        d['pragmas.' + key] = val
    dd.append(d)

df = pandas.DataFrame(dd)

ind = (df.data == 'gencode.v19.annotation.gtf.10k')# & (df.infer_gene_extent)
wd = df.ix[ind].sort(['write'])
rd = wd



fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(wd.write, color='k', linewidth=2)
ax2.plot(rd.read, color='k', linewidth=2)


def add_symbols(key, val, style, offset=0):
    wxs = []
    rxs = []
    label='%s=%s' % (key, val)
    _style = dict(
        linestyle='None',
        label=label,
        alpha=0.5)
    _style.update(style)

    for i, (_, row) in enumerate(wd.iterrows()):
        if row[key] == val:
            wxs.append(i)
    for i, (_, row) in enumerate(rd.iterrows()):
        if row[key] == val:
            rxs.append(i)

    (xmin, xmax, ymin, ymax) = ax1.axis()
    _offset = offset * ((ymax - ymin) / 10)
    line, = ax1.plot(wxs, [ymin+_offset for _ in wxs], **_style)
    (xmin, xmax, ymin, ymax) = ax2.axis()
    _offset = offset * ((ymax - ymin) / 10)
    line, ax2.plot(rxs, [ymin+_offset for _ in rxs], **_style)
    return line, label


for_legend = []
offset = 0
for (column, symbol) in [
    ('pragmas.journal_mode', 'o'),
    ('pragmas.locking_mode', '^'),
    ('pragmas.temp_store', 's'),
    ('pragmas.synchronous', '>'),
]:
    offset += 0.5
    for option in rd[column].unique():
        for_legend.append(add_symbols(column, option, dict(marker=symbol), offset=offset))


lines, labels = zip(*for_legend)
fig.legend(lines, labels, loc='center right', numpoints=1, prop=dict(size=8))

ax1.set_ylabel('write time (s)')
ax2.set_ylabel('read time (s)')
ax1.set_xlabel('rank')
ax2.set_xlabel('rank')
fig.tight_layout()
fig.subplots_adjust(right=0.7)
fig.savefig('runs.png')

fig = plt.figure(figsize=(8, 10))
cols = [
    'pragmas.journal_mode',
    'pragmas.temp_store',
    'pragmas.synchronous',
    'pragmas.locking_mode']
ncols = len(cols)
nrows = 2
i = 0
for col in cols:
    i += 1
    ax = fig.add_subplot(ncols, 2, i)
    sns.boxplot(
        rd['write'],
        groupby=rd[col],
        ax=ax)
    i += 1
    ax = fig.add_subplot(ncols, 2, i)
    sns.boxplot(
        rd['read'],
        groupby=rd[col],
        ax=ax)
fig.tight_layout()
fig.savefig('boxplots.png')
plt.show()
