#!/usr/bin/env python

from setuptools import setup, find_packages

install_requires = set()
with open( "requirements.txt" ) as f:
  for dep in f.read().split('\n'):
      if dep.strip() != '' and not dep.startswith('-e'):
          install_requires.add( dep )

setup(name='MODIS_Aggregation',
      version='0.0.1',
      description='Creates classes for MODIS aggregation',
      author='Jianwu Wang',
      zip_safe=False,
      author_email='jianwu@umbc.edu',
      url='https://github.com/big-data-lab-umbc/MODIS_Aggregation',
      packages=find_packages(),
      install_requires=list(install_requires),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Apache 2.0 License",
        "Operating System :: OS Independent",
    ],
)

