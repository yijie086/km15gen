from setuptools import setup
from Cython.Build import cythonize

setup(
    name='km15gen',
    ext_modules=cythonize("km15gen.pyx"),
    zip_safe=False,
)

