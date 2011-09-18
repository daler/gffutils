from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
        name='gffutils',
        cmdclass={'build_ext': build_ext},
        ext_modules=[Extension('gfffeature', sources=['gffutils/gfffeature.pyx'])]
    )
