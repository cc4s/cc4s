# -*- coding: utf-8 -*-

import sys
import glob
from setuptools import setup, find_packages

__version__ = "0.0.1"
__author__ = "Alejandro Gallo"
__email__ = "aamsgallo@gmail.com"
__license__ = "GPLv3"


setup(
    name='testis',
    version=__version__,
    maintainer='Alejandro Gallo',
    maintainer_email='aamsgallo@gmail.com',
    author=__author__,
    author_email=__email__,
    license=__license__,
    url='https://github.com/testis/testis',
    install_requires=[
        "PyYAML",
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'testis=testis:main',
        ],
    },
    platforms=['linux', 'osx', 'win32'],
)
