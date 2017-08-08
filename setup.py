
import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'gffutils', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')
requirements = open(os.path.join(os.path.dirname(__file__), 'requirements.txt')).readlines()
setup(
    name='gffutils',
    version=version,
    install_requires=requirements,
    packages=['gffutils', 'gffutils.scripts', 'gffutils.test',
              'gffutils.test.data'],
    scripts=['gffutils/scripts/gffutils-cli'],
    author='Ryan Dale',
    package_dir={'gffutils': 'gffutils'},
    package_data = {'gffutils': ['test/data/*']},
    description="Work with GFF and GTF files in a flexible "
    "database framework",
    author_email='dalerr@niddk.nih.gov',
    url='https://github.com/daler/gffutils',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
