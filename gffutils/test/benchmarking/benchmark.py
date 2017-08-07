import os
import gffutils
import time
import itertools
import pandas
import sqlite3
import json
import subprocess

if not os.path.exists('timings.db'):
    conn = sqlite3.connect('timings.db')
    c = conn.cursor()
    c.executescript(
        """
        CREATE TABLE timings (
        branch text,
        hash text,
        write float,
        read float,
        kwargs text,
        primary key (branch, hash, kwargs));
        """)
    conn.commit()

conn = sqlite3.connect('timings.db')

branches = subprocess.check_output(['git', 'branch'])
branch = [i for i in branches.splitlines(False) if i.startswith('*')]
branch = branch[0].replace('* ', '')
commit = subprocess.check_output(['git', 'describe']).strip()

class Benchmark(object):
    """
    Class for benchmarking writing and reading gffutils dbs.
    """

    def timeit(method):
        def timed(*args, **kw):
            self = args[0]
            ts = time.time()
            result = method(*args, **kw)
            te = time.time()
            self.results[method.__name__] = te-ts
            return result
        return timed

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.serialized_kwargs = json.dumps(kwargs)
        self.results = {}
        self.db = None

    def run(self):
        res = self.exists()
        if res:
            return res
        else:
            self.write()
            self.read()
            print(dict(self.results))
            self.insert()
        return

    def insert(self):
        d = dict(self.results)
        c = conn.cursor()
        c.execute(
            """
            INSERT OR IGNORE INTO timings VALUES (?, ?, ?, ?, ?)
            """,
            (
                branch,
                commit,
                self.results['write'],
                self.results['read'],
                self.serialized_kwargs
            )
        )
        conn.commit()

    def exists(self):
        c = conn.cursor()
        res = c.execute('select * from timings where branch = ? and hash = ? and kwargs = ?',
                        (branch, commit, self.serialized_kwargs))
        return list(res)

    @timeit
    def write(self):
        self.db = gffutils.create_db(**self.kwargs)


    def get_db(self):
        if self.db is None:
            raise ValueError("Run create() first")
        return self.db


    @timeit
    def read(self):
        db = self.get_db()
        for i in db.features_of_type('gene'):
            for j in db.children(i):
                pass
        for i in db.features_of_type('exon'):
            for j in db.parents(i):
                pass



base_kwargs = dict(
    data='gencode.v19.annotation.gtf.10k',
    dbfn=':memory:',
    force=True, infer_gene_extent=False)


pragmas = {
    'journal_mode': ['MEMORY', 'DELETE', 'WAL', 'OFF'],
    'synchronous': ['OFF', 'NORMAL'],
    'main.cache_size': [2000, 10000, 100000],
    'page_size': [1024, 4096, 8192],
    'temp_store': [0, 1, 2],
    'locking_mode': ['NORMAL', 'EXCLUSIVE'],
}

varNames = sorted(pragmas)
combinations = [dict(zip(varNames, prod)) for prod in itertools.product(*(pragmas[varName] for varName in varNames))]
for c in combinations:
    for k, v in c.items():
        if v is None:
            del c[k]

combinations = sorted([dict(i) for i in set([tuple(c.items()) for c in combinations])])

d = []

n = len(combinations)
k = 0
for pragma in combinations:
    k += 1
    kwargs = base_kwargs.copy()
    kwargs['pragmas'] = pragma
    kwargs['infer_gene_extent'] = True
    b = Benchmark(**kwargs)
    print('%s of %s' % (k, n))
    res = b.run()
conn.commit()

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

df_exp = pandas.DataFrame(dd)
