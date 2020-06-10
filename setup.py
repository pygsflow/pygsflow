#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script. """

from setuptools import setup, find_packages
import sys

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

if sys.version_info >= (3, 0):
    requirements = ["pandas",
                    "numpy",
                    "flopy >= 3.3.1",
                    "pyshp",
                    "pycrs",
                    "matplotlib"]
else:
    requirements = ["pandas",
                    "numpy",
                    "flopy == 3.2.13",
                    "pyshp",
                    "pycrs",
                    "matplotlib"]

setup_requirements = []

test_requirements = []

setup(
    author="Ayman Alzraiee, Joshua Larsen",
    author_email='aalzraiee@usgs.gov, jlarsen@usgs.gov',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Hydrology'
    ],
    description="""Python package to read, edit, write, and visualize GSFLOW models.
                   This package relies on FloPy for reading modflow input, but adds additional
                   fuctionality and flexibility through it's interface. It includes additional
                   support for the new MODFLOW-NWT AG package. The plotting library is also 
                   compatible with FloPy's PlotMapView and PlotCrossSection.""",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='gsflow',
    name='gsflow',
    packages=find_packages(include=['gsflow',
                                    'gsflow.prms',
                                    'gsflow.utils',
                                    'gsflow.modflow',
                                    'gsflow.modsim',
                                    'gsflow.output']),
    setup_requires=setup_requirements,
    test_suite='autotest',
    tests_require=test_requirements,
    url='https://github.com/aymanalz/gsflow',
    version='0.2.0',
    zip_safe=False,
)
