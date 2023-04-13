from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import glob



extensions = [
    Extension("processing_module", sources=["processing_module.pyx"], include_dirs=[numpy.get_include()], extra_compile_args=["-O3"], language="c++")
]

setup(
    name='km15gen_c',
    ext_modules=cythonize(["km15gen.pyx"] , language_level = 3, annotate=True),
    zip_safe=False
)