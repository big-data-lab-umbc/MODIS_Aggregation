#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages


install_requires = set()
with open("requirements.txt") as f:
    for dep in f.read().split('\n'):
        if dep.strip() != '' and not dep.startswith('-e'):
            install_requires.add(dep)

setup(
    author="Jianwu Wang",
    author_email='jianwu@umbc.edu',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Modis",
    entry_points={
        'console_scripts': [
            'MODIS=MODIS.cli:main',
        ],
    },
    install_requires=list(install_requires),
    license="Apache Software License 2.0",
    long_description="MODIS AGGREGATION",
    include_package_data=True,
    keywords='MODIS_Aggregation',
    name='MODIS_Aggregation',
    packages=find_packages(),
    test_suite='tests',
    url='https://github.com/big-data-lab-umbc/MODIS_Aggregation',
    version='0.1.0',
    zip_safe=False,
)
