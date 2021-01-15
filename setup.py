import os
import platform
import subprocess
import sys

from setuptools import setup
from setuptools.extension import Extension

import numpy as np

name = os.path.basename(os.path.abspath('.'))
print(name)

with open('VERSION', 'r') as f:
    __version__ = f.read().strip()

def query_config(config, query):
    ret = subprocess.run('{} {}'.format(config, query).split(), stdout=subprocess.PIPE)
    return str(ret.stdout, 'ascii').strip()

fj_cxxflags = query_config('fastjet-config', '--cxxflags')
cxxflags = ['-fopenmp', '-std=c++14', '-DSWIG_TYPE_TABLE=fastjet', fj_cxxflags]
ldflags = ['-fopenmp']
libs = [name]
if platform.system() == 'Darwin':
    cxxflags.insert(0, '-Xpreprocessor')
    del ldflags[0]
    libs.append('omp')

fj_lib_flags = query_config('fastjet-config', '--libs').split()
fj_libdirs = [x[2:] for x in fj_lib_flags if x.startswith('-L')]
libs += [x[2:] for x in fj_lib_flags if x.startswith('-l')]

if sys.argv[1] == 'swig':

    fj_prefix = query_config('fastjet-config', '--prefix')
    fj_swig_interface = os.path.join(fj_prefix, 'share', 'fastjet', 'pyinterface', 'fastjet.i')

    opts = '-fastproxy {} -DSWIG_NUMPY -DFASTJET_SWIG_INTERFACE={}'.format(fj_cxxflags, fj_swig_interface)
    if len(sys.argv) > 2:
        opts += ' ' + ' '.join(sys.argv[2:])

    command = 'swig -python -c++ {opts} -IWasserstein/wasserstein -o Py{name}.cc {name}.i'.format(opts=opts, name=name)
    print(command)
    subprocess.run(command.split())

else:

    with open('README', 'r') as f:
        readme = f.read()

    module = Extension('_' + name,
                   sources=['Py{}.cc'.format(name)],
                   language='c++',
                   include_dirs=[np.get_include(), 'Wasserstein/wasserstein'],
                   library_dirs=fj_libdirs,
                   libraries=libs,
                   extra_compile_args=cxxflags,
                   extra_link_args=ldflags,
                   #undef_macros=['NDEBUG']
                  )

    setup(
        name=name,
        version=__version__,
        author='Patrick T. Komiske III',
        description='{} FastJet Contrib'.format(name),
        long_description=readme,
        url='https://fastjet.hepforge.org/contrib/',
        py_modules=[name],
        ext_modules=[module]
    )
