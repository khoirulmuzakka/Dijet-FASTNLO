#! /usr/bin/env python

from distutils.core import setup, Extension

longdesc = """This is a simple SWIG wrapper on C++ fastNLO Reader interface."""

## Extension definition
import os
wrapsrc = './fastnlo_wrap.cc'
incdir_src = os.path.abspath('../include')
incdir_build = os.path.abspath('../include')
libdir = os.path.abspath('../fastnlotoolkit')

incdir_lhapdf = '/home/aem/uni/sw/lhapdf/include'
libdir_lhapdf = '/home/aem/uni/sw/lhapdf/lib'
libdir_fnlo = '/home/aem/uni/sw/fnlo_toolkit/lib'

cxxargs = '-g -O2'.split()
ldargs = ''.split()
ext = Extension('_fastnlo',
                [wrapsrc],
                include_dirs=[incdir_src, incdir_build,incdir_lhapdf],
                library_dirs=[libdir,libdir_lhapdf, os.path.join(libdir,'.libs')],
                runtime_library_dirs=[libdir_lhapdf,libdir_fnlo],
                extra_compile_args = cxxargs,
                extra_link_args = ldargs,
                libraries=['fastnlotoolkit','LHAPDF'],
                )

## Setup definition
setup(name = 'fastnlo',
      version = '2.1.0',
      #include_package_data = True,
      ext_modules=[ext],
      py_modules = ['fastnlo'],
      author = ['fastNLO'],
      url = 'http://projects.hepforge.org/fastnlo/',
      description = 'A SWIG wrapper on C++ fastNLO interface.',
      long_description = longdesc,
      )

