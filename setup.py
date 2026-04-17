# -*- coding: utf-8 -*-
"""setup.py: setuptools control for backward compatibility.
See:
https://packaging.python.org/en/latest/distributing.html
"""

from setuptools import setup, find_packages
import os
import re

# Extract version directly from the source file
def read_version():
    with open('tfbs_footprinter3/tfbs_footprinter3.py', encoding='utf-8') as f:
        return re.search(r'^__version__\s*=\s*"(.*)"', f.read(), re.M).group(1)

# Read the README for the long description
def read_readme():
    with open('README.md', encoding='utf-8') as f:
        return f.read()

setup(
    name='TFBS_footprinting3',
    version=read_version(),
    description='Tool for identifying conserved TFBSs in vertebrate species.',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/thirtysix/TFBS_footprinting3',
    author='Harlan Barker',
    author_email='harlan.barker@tuni.fi',
    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],

    keywords='bioinformatics analysis',
    packages=find_packages(include=['tfbs_footprinter3']),
    install_requires=[
        'httplib2',
        'numpy',
        'matplotlib',
        'biopython',
        'msgpack',
        'wget',
        'pandas'
    ],
    extras_require={},
    entry_points={
        'console_scripts': [
            'tfbs_footprinter3 = tfbs_footprinter3.tfbs_footprinter3:main'
        ]
    },
    include_package_data=True,
    package_dir={'tfbs_footprinter3': 'tfbs_footprinter3'},
    package_data={'tfbs_footprinter3': ['data/*.json']},
)
