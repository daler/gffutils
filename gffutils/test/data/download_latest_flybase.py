import os
import sys
import time
import argparse
import ftplib
import re
import gffutils

usage = """


"""
CHROMS = ['2L', '2LHet', '2R', '2RHet', '3L', '3LHet', '3R', '3RHet', '4', 'U',
'Uextra', 'X', 'Xhet', 'YHet', 'all-no-analysis', 'all',
'dmel_mitochondrion_genome']

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument('--dest', default=".", help="Destination directory; "
        "FlyBase GFF will be downloaded here and database will be created")
ap.add_argument('--list', action="store_true", help="Just list the available"
        "chromosomes to download")
ap.add_argument('--chrom', default='all-no-analysis', help="chromosome"
        " or version to download, one of %s" % CHROMS)
ap.add_argument('--make_db', action='store_true',
        help="Create a gffutils database from the downloaded file")
args = ap.parse_args()


URL = "ftp.flybase.net"
PATH = "genomes/dmel/current/gff/"

ftp = ftplib.FTP(URL)
sys.stderr.write(ftp.login())

dirlist = ftp.nlst(PATH)
sys.stderr.write('\n')
for i in dirlist:
    size = ftp.size(i)
    if ('-' + args.chrom + '-') in i:
        fn = i
        flag = 'X'
    else:
        flag = " "

    sys.stderr.write("[ %s ] [%.1f MB] %s\n" % (flag, (size / 1e6), i))

size = ftp.size(fn)

class Progress(object):
    def __init__(self, fileobj, size, template="%s%%"):
        """
        blocksize is what retrbinary uses; default is 8192
        """
        self.fileobj = fileobj
        self.size = size
        self.total = 0.0
        self.last_perc = 0
        self.template = template

    def __call__(self, data):
        if size is not None:
            self.total += len(data)
            perc = self._perc()
            if perc != self.last_perc:
                sys.stderr.write('\r' + (self.template % perc))
                sys.stderr.flush()
            self.last_perc = perc

        self.fileobj.write(data)

    def _perc(self):
        return int(self.total / self.size * 100)


if not args.list:
    t0 = time.time()
    f = open(os.path.join(args.dest, os.path.basename(fn)), 'w')
    P = Progress(f, size)
    ftp.retrbinary("RETR " + fn, P)
    f.close()
    t1 = time.time()
    sys.stderr.write(' (%.2fMB done in %.2fs)\n' % ((size / 1e6), (t1 - t0)))

    if args.make_db:
        gffutils.create_db(f.name, f.name + '.db', verbose=True, force=True)
