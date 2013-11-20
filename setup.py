#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext
import numpy
from numpy.distutils.system_info import get_info


mpfit_sources = [
    'mpyfit/mpyfit.c',
    'mpyfit/cmpfit/mpfit.c',
]
# Avoid some numpy warnings, dependent on the numpy version
npy_api_version = "NPY_{0}_{1}_API_VERSION".format(
    *(numpy.__version__.split('.')[:2]))

setup(
    name='mpyfit',
    version='0.9.0',
    license='BSD 2 part license',
    maintainer='Evert Rol',
    maintainer_email='e.rol@sron.nl',
    packages=['mpyfit'],
    url='https://github.com/evertrol/mpyfit',
    classifiers=[
        'Intentended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programing Language :: Python',
        'Licence :: OSI Approved :: BSD License',
    ],
    ext_modules=[
        Extension(
            'mpyfit.mpfit',
            sources=mpfit_sources,
            include_dirs=['mpyfit/cmpfit', numpy.get_include()],
            extra_compile_args=['-std=c99 '
                                '-DNPY_NO_DEPRECATED_API={0}'.format(
                    npy_api_version)]
        ),
    ]
)
