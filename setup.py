#! /usr/bin/env python
import sys
import os
from distutils.core import setup, Extension

srcFileList= ["spilltree.i"]
moduleName= "spilltree2D" # may be reset from setup.cfg below

# Scan through the setup.cfg file to get the module name,
# because we have to make the shared object name match.
f= open("setup.cfg","r")
lines= f.readlines()
f.close()
for line in lines:
    words= line.strip().split()
    for i in xrange(len(words)):
        if words[i] == "-module":
            moduleName= words[i+1]

setup(name=moduleName,
      version='1.0.7',
      ext_modules=[Extension('_'+moduleName,srcFileList,
                             include_dirs=['.'],
                             define_macros=[('NDEBUG','1')],
                             #extra_compile_args=['-g','-O0'],
                             libraries=[],
                             swig_opts=['-c++'],
                             ),
                   ],
      py_modules=[moduleName],
      description='An implementation of Spilltrees, for finding k nearest neighbors',
      author='Joel Welling',
      author_email='welling@psc.edu'
      )
