from distutils.core import setup
from Cython.Build import cythonize

setup(name = 'mfd',
      version = '0.1',
      url = 'https://github.com/UP-RS-ESP/mfdrouting/cython',
      author = 'Aljoscha Rheinwalt',
      author_email = 'aljoscha.rheinwalt@uni-potsdam.de',
      ext_modules = cythonize('mfd.pyx')
)
