"""
Benchmarking script, trying different methods of maintaining order when
populating a dictionary.
"""
import time
import sys
from gffutils.helpers import DefaultOrderedDict
from gffutils.feature import Attributes
from collections import defaultdict
from collections import OrderedDict
import random


nlines = 200000
nkeys = 5
vals_per_key = 3


seen = set()
order = []

def attributesupdate(d, k, v, seen, order):
    """
    Append value `v` to list at key `k` in dict `d`
    """
    val = d.setdefault(k, [])
    val.append(v)
    return seen, order

def dictupdate(d, k, v, seen, order):
    """
    Append value `v` to list at key `k` in dict `d`.
    """
    val = d.setdefault(k, [])
    val.append(v)
    return seen, order


def defaultupdate(d, k, v, seen, order):
    """
    Append value `v` to list at key `k` in dict `d`.
    """
    d[k].append(v)
    return seen, order

def naiveorderupdate(d, k, v, seen, order):
    defaultupdate(d, k, v, seen, order)
    if k not in seen:
        seen.update([k])
        order.append(k)
    return seen, order

lookup = {
    'dict': (lambda: dict(), dictupdate),
    'defaultordered': (lambda: DefaultOrderedDict(list), defaultupdate),
    'default': (lambda: defaultdict(list), defaultupdate),
    'ordered': (lambda: OrderedDict(), dictupdate),
    'naiveordered': (lambda: defaultdict(list), naiveorderupdate),
    'attributes': (lambda: Attributes(), attributesupdate),
}

try:
    kinds = sys.argv[1:]
except IndexError:
    sys.stderr.write("usage: %s [%s]\n" % (sys.argv[0], '|'.join(lookup.keys())))
    sys.exit(1)

for kind in kinds:
    func, update = lookup[kind]

    t0 = time.time()
    for line in range(nlines):
        d = func()
        seen = set()
        order = []
        for i in range(nkeys):
            for j in range(vals_per_key):
                seen, order = update(d, i, j, seen, order)
    print "%s: %.3fs" % (kind, time.time() - t0)

# check at least the last created dict to make sure it's correct
#assert dict(d) == {0: [0, 1, 2], 1: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2], 4: [0, 1, 2]}
