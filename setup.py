from setuptools import setup
from distutils.extension import Extension
try:
    from Cython.Distutils import build_ext
except ImportError:
    print "-------------------------------------------------"
    print "Please install Cython: (try 'pip install cython')"
    print "-------------------------------------------------"

setup(
        name='gffutils',
        cmdclass={'build_ext': build_ext},
        install_requires=['cython'],
        ext_modules=[Extension('gfffeature', sources=['gffutils/gfffeature.pyx'])],
        packages=['gffutils'],
        author='Ryan Dale',
        package_dir={'gffutils': 'gffutils'},
        description="Work with GFF and GTF files in a flexible database framework",
        author_email='dalerr@niddk.nih.gov',
        url='none',
    )
