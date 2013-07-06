import gffutils
dialect = gffutils.constants.dialect
fn = gffutils.example_filename('dmel-all-no-analysis-r5.49_50k_lines.gff')
#for line in open(fn):
#    gffutils.parser._split_keyvals(line.split('\t')[-1], dialect=dialect)
db = gffutils.create_db(fn, ':memory:', merge_strategy='merge')
