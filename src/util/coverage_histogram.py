#!/usr/bin/env python

# expected input:
# - genome<tab>depth<tab>bases
# e.g. 
# bedtools genomecov -ibam {mapped_bam} | python coverage_histogram plot.html > stats.txt

import collections
import sys

import plotly
import sys

from plotly.graph_objs import Scatter, Layout, Bar

c = 0
t = 0
cex0 = 0
tex0 = 0
h = collections.defaultdict(int)

for count, l in enumerate(sys.stdin):
        f = l.split('\t')
        if f[0] != 'genome':
                continue
        depth = int(f[1])
        bases = int(f[2])
        c += depth * bases
        t += bases
        if depth > 0:
                cex0 += depth * bases
                tex0 += bases
        h[depth] += bases
        if count % 10000000 == 0:
            sys.stderr.write('coverage_histogram: processed {} lines\n'.format(count))

sys.stdout.write('coverage_total: {}\ncoverage_bases: {}\ncoverage_mean: {}\ncoverage_total_no_zero: {}\ncoverage_mean_no_zero: {}'.format(c, t, 1.0 * c/t, cex0, 1.0 * cex0/tex0))

x = [x for x in range(1,100)]

plotly.offline.plot({
    "data": [Scatter(x=x, y=[h[v] for v in x])],
    "layout": Layout(title="Coverage", xaxis=dict(title='Coverage'), yaxis=dict(title='Count')),
    },
    filename=sys.argv[1],
    auto_open=False)
