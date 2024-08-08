#!/usr/bin/env python

from distutils.core import setup
from setuptools import setup, find_packages



setup(name='TACHELES',
      version='0.1',
      description='A transit model for exoplanets and stars with extended outer layers',
      author='Volker Perdelwitz',
      author_email='volker.perdelwitz@weizmann.ac.il',
      url='https://www.weizmann.ac.il/EPS/',
      # project_urls={'Source code': '"'https://github.com/pypa/sampleproject/issues"'}'
      keywords='TACHELES, transit, exoplanet, chromosphere, activity',
      package_dir = {'': 'src'},
      packages=find_packages(where='src'),
      # install_requires=['PyAstronomy','astropy', 'math', 'numba', 'numpy', 'photutils', 'scipy', 'batman', 'csv', 'matplotlib', 'os', 'time']
     )

