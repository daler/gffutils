import ez_setup
ez_setup.use_setuptools()

import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'gffutils', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')

setup(
    name='gffutils',
    version=version,
    install_requires=['argh', 'argcomplete', 'simplejson'],
    packages=['gffutils', 'gffutils.scripts', 'gffutils.test',
              'gffutils.test.data'],
    scripts=['gffutils/scripts/gffutils-cli'],
    author='Ryan Dale',
    package_dir={'gffutils': 'gffutils'},
    package_data = {'gffutils': ['test/data/*']},
    description="Work with GFF and GTF files in a flexible "
    "database framework",
    author_email='dalerr@niddk.nih.gov',
    url='none',
)
