import os
from setuptools import setup


def read(fname):
    data = ''
    with open(os.path.join(os.path.dirname(__file__), fname)) as infile:
        data = infile.read()
    return data

setup(
    name='impy',
    version='0.2',
    author="Shaman Narayanasamy, Yohan Jarosz",
    author_email="shaman.narayanasamy@uni.lu, yohan.jarosz@uni.lu",
    description=("Integrated Metaomic Pipeline command-line utility."),
    long_description=read('README.rst'),
    keywords="metagenomics metatranscriptomics pipeline integrated ecosystems biology",
    url="http://r3lab.uni.lu/web/imp",
    license = "MIT",
    py_modules=['impy'],
    install_requires=[
    ],
    entry_points='''
        [console_scripts]
        impy=impy:cli
    ''',
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License"
    ],
)
