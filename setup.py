#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
try:
    from setuptools import setup
    from setuptools.extension import Extension
    from setuptools.dist import Distribution
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    from distutils.dist import Distribution


ext_modules = []
if (any('--' + opt in sys.argv for opt in Distribution.display_option_names +
       ['help-commands', 'help']) or sys.argv[1] in ('egg_info', 'clean')):
    pass
else:
    import numpy
    mpfit_sources = [
        'mpyfit/mpyfit.c',
        'mpyfit/cmpfit/mpfit.c',
    ]
    # Avoid some numpy warnings, dependent on the numpy version
    npy_api_version = "NPY_{0}_{1}_API_VERSION".format(
        *(numpy.__version__.split('.')[:2]))
    ext_modules.append(
        Extension(
            'mpyfit.mpfit',
            sources=mpfit_sources,
            include_dirs=['mpyfit/cmpfit', numpy.get_include()],
            extra_compile_args=['-std=c99',
                                '-Wno-declaration-after-statement',
                                '-DNPY_NO_DEPRECATED_API={0}'.format(
                                    npy_api_version)]
        ),
    )


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
    requires=['numpy (>=1.6)'],
    ext_modules=ext_modules
)
