#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import imp
import json
import setuptools.command.develop
from setuptools import setup, find_packages
from importlib.util import module_from_spec, spec_from_file_location


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# Loads version.py module without importing the whole package.
def get_version_and_cmdclass(package_path):
    spec = spec_from_file_location('version',
                                   os.path.join(package_path, '_version.py'))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass


version, cmdclass = get_version_and_cmdclass('semicon')


def import_submodule(path, name):
    spec = spec_from_file_location(name, os.path.join(path, name + '.py'))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_cache(dir):
        print('building model cache')
        sys.path.append('semicon')
        from kp_models import explicit_foreman, explicit_zeeman
        sys.path.pop()
        data = {
            'foreman': str(explicit_foreman.foreman),
            'zeeman': str(explicit_zeeman.zeeman),
        }

        cache_file = os.path.join(dir, 'semicon', 'model_cache.json')
        with open(cache_file, 'w') as f:
            json.dump(data, f)


# Build model cache from 'kp_models' package
class build_py(cmdclass['build_py']):

    def run(self):
        # make sure we run the miniver stuff
        super().run()
        build_cache(self.build_lib)


class develop(setuptools.command.develop.develop):

    def run(self):
        super().run()
        data = build_cache('.')


cmdclass['build_py'] = build_py
cmdclass['develop'] = develop

classifiers = """\
    Development Status :: 3 - Alpha
    Intended Audience :: Science/Research
    Intended Audience :: Developers
    Programming Language :: Python :: 3 :: Only
    Topic :: Software Development
    Topic :: Scientific/Engineering
    Operating System :: POSIX
    Operating System :: Unix"""

setup(
    name="semicon",
    version=version,

    author='R.J. Skolasinski',
    author_email='r.j.skolasinski@gmail.com',
    description=("Package for simulating quantum mechanical k·p Hamiltonians"),
    license="BSD",

    long_description=read("README.md"),
    platforms=["Unix", "Linux"],
    url="https://gitlab.kwant-project.org/semicon/semicon",


    packages=find_packages('.'),
    package_data={'semicon': ['databank/*.yml']},

    setup_requires=['sympy >= 0.7.6'],
    install_requires=['pyyaml', 'scipy >= 0.17', 'kwant >= 1.3',
                      'sympy >= 1.1.1', 'pandas >= 0.19.2'],
    classifiers=[c.strip() for c in classifiers.split('\n')],
    cmdclass=cmdclass,
)
