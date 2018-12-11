#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import imp
import json
from setuptools import setup, find_packages


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# Building cache

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
cache_fname = os.path.join(BASE_DIR, 'semicon', 'model_cache.json')
csv_fnames = os.path.join(BASE_DIR, 'semicon', 'databank', '*.csv')


def build_cache():
    print("building models' cache")
    from kp_models import explicit_foreman, explicit_zeeman

    data = {
        'foreman': str(explicit_foreman.foreman),
        'zeeman': str(explicit_zeeman.zeeman),
    }
    with open(cache_fname, 'w') as f:
        json.dump(data, f)


# Standard python build
print("path", BASE_DIR)

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
    version="0.0.0",

    author='R.J. Skolasinski',
    author_email='r.j.skolasinski@gmail.com',
    description=("Python package for doing k·p simulation"),
    license="BSD",

    long_description=read("README.md"),
    platforms=["Unix", "Linux"],
    url="https://gitlab.kwant-project.org/semicon/semicon",


    packages=find_packages('.'),
    package_data={'': [cache_fname, csv_fnames]},

    setup_requires=['sympy >= 0.7.6'],
    install_requires=['kwant >= 1.3', 'sympy >= 0.7.6', 'pandas >= 0.19.2'],
    classifiers=[c.strip() for c in classifiers.split('\n')]
)


build_cache()
