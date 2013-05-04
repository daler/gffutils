import ez_setup
ez_setup.use_setuptools()

import os
import sys
from setuptools import setup
from distutils.extension import Extension

##
## Why is Cython required for non-developers?  Can distribute C code...
##
try:
    from Cython.Distutils import build_ext

##
## Very few people can actually run this successfully:
##
##   easy_install -U cython
##
## virtually on all systems, setuptools is installed
## and not distribute, so it will just trigger a
## '-U option not recognized' error.
##
except ImportError:
    sys.stderr.write("""
==================================================

Please install Cython (http://cython.org/),
which is required to build pybedtools. Usually
you can do:

    pip install -U cython

or

    easy_install -U cython

==================================================
    """)
    sys.exit(1)

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'gffutils', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')

setup(
    name='gffutils',
    version=version,
    cmdclass={'build_ext': build_ext},
    install_requires=['cython', 'argh', 'argcomplete'],
    ext_modules=[
        Extension(
            'gffutils.gfffeature',
            sources=['gffutils/gfffeature.pyx'])
    ],
    packages=['gffutils', 'gffutils.scripts'],
    scripts=['gffutils/scripts/gffutils-cli'],
    author='Ryan Dale',
    package_dir={'gffutils': 'gffutils'},
    description="Work with GFF and GTF files in a flexible "
    "database framework",
    author_email='dalerr@niddk.nih.gov',
    url='none',
)
