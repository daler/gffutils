import os
import sys
from setuptools import setup
from distutils.extension import Extension

# optional cython
try:
  from Cython.Distutils import build_ext
except ImportError:
  from distutils.command import build_ext as _build_ext
  class build_ext(_build_ext.build_ext):

      description = "change pyx files to corresponding .c/.cpp (fallback when cython is not installed)"

      def build_extensions(self):
          # First, sanity-check the 'extensions' list
          self.check_extensions_list(self.extensions)

          for extension in self.extensions:
              iscpp = extension.language and extension.language.lower() == 'c++'
              target_ext = '.cpp' if iscpp else '.c'

              patchedsrc = []
              for source in extension.sources:
                (root, ext) = os.path.splitext(source)
                if ext == '.pyx':
                  patchedsrc.append(root + target_ext)
                else:
                  patchedsrc.append(source)

              extension.sources = patchedsrc
              self.build_extension(extension)


if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__


version_py = os.path.join(os.path.dirname(__file__), 'gffutils', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')

setup(
        name='gffutils',
        version=version,
        cmdclass={'build_ext': build_ext},
        install_requires=['cython'],
        ext_modules=[Extension('gffutils.gfffeature', sources=['gffutils/gfffeature.pyx'])],
        packages=['gffutils'],
        author='Ryan Dale',
        package_dir={'gffutils': 'gffutils'},
        description="Work with GFF and GTF files in a flexible database framework",
        author_email='dalerr@niddk.nih.gov',
        url='none',
    )
