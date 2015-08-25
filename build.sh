#!/usr/bin/env bash
echo "Compiling..."
python -c "from distutils.core import setup;from Cython.Build import cythonize;setup(ext_modules = cythonize('*.pyx'))" build_ext --inplace
echo "Done."