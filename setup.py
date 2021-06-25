#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script. """

from setuptools import setup, find_packages
import sys
import os


with open('README.md') as readme_file:
    readme = readme_file.read()

if sys.version_info >= (3, 0):
    requirements = ["pandas",
                    "numpy",
                    "flopy >= 3.3.1",
                    "pyshp",
                    "pycrs",
                    "matplotlib"]
else:
    raise EnvironmentError("pyGSFLOW is only supported with python 3")

setup_requirements = []

test_requirements = []

try:
    import pypandoc

    fpth = os.path.join(".", "RELEASE.md")
    long_description = pypandoc.convert_file(fpth, 'rst')
except ImportError:
    long_description = ""


setup(
    author="Ayman Alzraiee, Joshua Larsen, Rich Niswonger",
    author_email='aalzraiee@usgs.gov, jlarsen@usgs.gov, rniswon@usgs.gov',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Hydrology',
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.4",
    description="pyGSFLOW is a python package to create, run, and " +
                "post-process GSFLOW-based models",
    install_requires=requirements,
    license="MIT license",
    long_description=long_description,
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
    url='https://github.com/pygsflow/pygsflow',
    version='1.0.1',
    zip_safe=False,
)
