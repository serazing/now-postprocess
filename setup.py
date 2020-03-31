#!/usr/bin/env python

from os.path import exists
from setuptools import setup

DISTNAME = 'now-postprocess'
PACKAGES = ['nowpp']
TESTS = [p + '.tests' for p in PACKAGES]
INSTALL_REQUIRES = ['xarray', 'numpy >= 1.11',
                    'matplotlib>=1.5',
                    'pandas',]
TESTS_REQUIRE = ['pytest >= 2.7.1']
URL = ''
AUTHOR = 'Guillaume Serazin'
AUTHOR_EMAIL = 'guillaume.serazin@unsw.edu.au'
LICENSE = 'MIT'
DESCRIPTION = 'Collection of funtions to process NOW outputs'
LONG_DESCRIPTION = (open('README.md').read() if exists('README.md') else '')
VERSION = 0.1

SCRIPTS  = ['bin/calculate_mslp_arguments.ncl',
            'bin/runncl.sh',
            'bin/nowpp-mslp',
            'bin/nowpp-get']

#DATA_FILES=[('config', ['cfg/config.cfg'])]

setup(name=DISTNAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      scripts=SCRIPTS,
      url=URL,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=PACKAGES + TESTS,
      install_requires=INSTALL_REQUIRES,
      zip_safe=False)
