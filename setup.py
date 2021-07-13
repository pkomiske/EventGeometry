# EventGeometry Package
#
#  Questions/comments? pkomiske@mit.edu
#
#  Copyright (c) 2019-2021
#  Patrick T. Komiske III, Eric M. Metodiev, Jesse Thaler
#
#----------------------------------------------------------------------
# This file is part of FastJet contrib.
#
# It is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# It is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------

from __future__ import print_function

import os
import platform
import re
import subprocess
import sys

################################################################################

# Package name, with capitalization
name = 'EventGeometry'

# Python package name, lower case by convention
lname = name.lower()

# some extra options that swig needs for this package
extra_swig_opts = '-keyword -w325,402,509,511 -IWasserstein -DSWIG_NUMPY'

# using PyFJCore or not
use_pyfjcore = True

################################################################################

# function to query a config binary and get the result
fastjet_config = os.environ.get('FASTJET_CONFIG', 'fastjet-config')
def query_config(query):
    if not use_pyfjcore:
        return subprocess.check_output([fastjet_config, query]).decode('utf-8').strip()
    return ''

# get fastjet info
fj_prefix = query_config('--prefix')
fj_cxxflags = query_config('--cxxflags')
fj_ldflags = query_config('--libs')

# get contrib README
with open('README.md', 'r') as f:
    readme = f.read()

# get contrib version
with open('VERSION', 'r') as f:
    __version__ = f.read().strip()

HELP_MESSAGE = """{name} FastJet Contrib Python Package

Usage: python3 setup.py [COMMAND] [OPTIONS]

Valid options for COMMAND include:
  help - Show this message
  swig - Run SWIG to generate new {lname}.py and Py{name}.cc files; OPTIONS are passed to SWIG
  build_ext - Build the Python extension
  install - Install the Python extension to a standard location
  clean - Remove Python build directories

OPTIONS are passed along to setuptools commands such as build_ext, install, and clean.

Relevant environment variables include:
  FASTJET_CONFIG - Path to fastjet-config binary [defaults to looking for fastjet-config in PATH]
  CXXFLAGS - Compiler flags passed along when building the Python extension module
"""

def show_help():
    print(HELP_MESSAGE.format(name=name, lname=lname))

def run_swig():

    contrib = {'docstring': '{}'.format(readme.replace('"', r'\"')),
               'version': "{}".format(__version__)}

    interface_file = '{}.i'.format(lname)
    template_file = '{}.i.template'.format(lname)
    print('Constructing SWIG interface file {} from {}'.format(interface_file, template_file))

    # read interface template and write interface file
    with open(template_file, 'r') as f_template, open(interface_file, 'w') as f_interface:
        temp = f_template.read().replace(r'\{', r'\<<').replace(r'\}', r'\>>')
        temp = temp.replace('{', '{{').replace('}', '}}').replace(r'\<<', '{').replace(r'\>>', '}')
        f_interface.write(temp.format(**contrib))

    # form swig options
    opts = '-fastproxy {}'.format(fj_cxxflags)
    if use_pyfjcore:
        opts += ' -DEVENTGEOMETRY_USE_PYFJCORE -IPyFJCore/pyfjcore'
    else:
        opts += ' -DFASTJET_PREFIX=' + fj_prefix

    # handle extra options for swig
    sys.argv += extra_swig_opts.split()
    if len(sys.argv) > 2:
        opts += ' ' + ' '.join(sys.argv[2:])

    command = 'swig -python -c++ {} -o Py{}.cc {}.i'.format(opts, name, lname)
    print(command)
    subprocess.run(command.split())

def run_setup():

    # get cxxflags from environment, add fastjet cxxflags, and SWIG type table info
    cxxflags = os.environ.get('CXXFLAGS', '').split() + fj_cxxflags.split() + ['-g0', '-fopenmp']
    libs, ldflags = [], []

    # handle multithreading with OpenMP
    if platform.system() == 'Darwin':
        cxxflags.insert(-1, '-Xpreprocessor')
        libs.append('omp')
    else:
        ldflags.append('-fopenmp')

    # determine library paths and names for Python
    library_dirs = []
    for x in fj_ldflags.split():
        if x.startswith('-L'):
            library_dirs.append(x[2:])
        elif x.startswith('-l'):
            libs.append(x[2:])
        else:
            ldflags.append(x)

    from setuptools import setup
    from setuptools.extension import Extension

    import numpy as np

    sources = ['Py{}.cc'.format(name)]
    include_dirs = [np.get_include(), 'Wasserstein']
    macros = [('SWIG_TYPE_TABLE', 'fastjet')]
    if use_pyfjcore:
        sources.append('PyFJCore/pyfjcore/fjcore.cc')
        cxxflags.append('-std=c++14')
        include_dirs.append('PyFJCore')
        macros.append(('EVENTGEOMETRY_USE_PYFJCORE', None))

    module = Extension('{0}._{0}'.format(lname),
                       sources=sources,
                       language='c++',
                       include_dirs=include_dirs,
                       library_dirs=library_dirs,
                       libraries=libs,
                       extra_compile_args=cxxflags,
                       extra_link_args=ldflags,
                       define_macros=macros
                      )

    setup(
        version=__version__,
        ext_modules=[module],
    )

def main():
    commands = {'help': show_help, 'swig': run_swig}
    commands.get(sys.argv[1], run_setup)()

if __name__ == '__main__':
    main()
