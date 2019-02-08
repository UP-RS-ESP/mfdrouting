from distutils.core import setup, Extension
import numpy

mod = Extension('mfdrouting',
    include_dirs = [numpy.get_include()],
    sources = ['mfdrouting.c'],
)

setup (name = 'mfdrouting',
    author = 'Aljoscha Rheinwalt',
    author_email = 'aljoscha.rheinwalt@uni-potsdam.de',
    ext_modules = [mod]
)
