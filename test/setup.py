# -*- coding: utf-8 -*-

import sys
import glob
from setuptools import setup, find_packages
import testis


setup(
    name='testis',
    version=testis.__version__,
    maintainer='Alejandro Gallo',
    maintainer_email='aamsgallo@gmail.com',
    author=testis.__author__,
    author_email=testis.__email__,
    license=testis.__license__,
    url='https://github.com/testis/testis',
    install_requires=[
        "PyYAML>=3.12",
        "tqdm>=4.1",
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'testis=testis:main',
        ],
    },
    platforms=['linux', 'osx', 'win32'],
)
