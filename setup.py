from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy



extensions = [
    Extension("processing_module", sources=["processing_module.pyx"], include_dirs=[numpy.get_include()], extra_compile_args=["-O3"], language="c++")
]


setup(
    name='km15gen',
    ext_modules=cythonize("km15gen.pyx"),
    zip_safe=False,
)

