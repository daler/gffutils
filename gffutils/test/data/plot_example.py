"""
Example to plot random genes
"""

import gffutils
from gffutils.contrib.plotting import Gene
from pylab import *

G = gffutils.FeatureDB('dmel-all-no-analysis-r5.43.gff.db')
gene = G.random_feature('gene')

gene_collection = Gene(
        G, 
        gene,
        utrs=['three_prime_UTR', 'five_prime_UTR'],
        color="0.5",
        edgecolor="None")

fig = figure()
ax = fig.add_subplot(111)
gene_collection.add_to_ax(ax)
ax.axis('tight')
show()
