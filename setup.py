import os
from setuptools import setup


def read(fname):
    data = ''
    with open(os.path.join(os.path.dirname(__file__), fname)) as infile:
        data = infile.read()
    return data

setup(
    name='impy',
    version=read('VERSION'),
    author="Shaman Narayanasamy, Yohan Jarosz",
    author_email="shaman.narayanasamy@uni.lu, yohan.jarosz@uni.lu",
    description=("Integrated Metaomic Pipeline command-line utility."),
    long_description=read('README.md'),
    keywords="metagenomics meatatranscriptomics pipeline integrated",
    url="http://r3lab.uni.lu/web/imp",
    license = "<>",
    py_modules=['impy'],
    install_requires=[
        'Click',
        'Path.py'
    ],
    entry_points='''
        [console_scripts]
        impy=impy:cli
    ''',
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: <> License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License"
    ],
)

setup(
    name = "an_example_pypi_project",
    version = "0.0.4",
    author = "Andrew Carter",
    author_email = "andrewjcarter@gmail.com",
    description = ("An demonstration of how to create, document, and publish "
                                   "to the cheese shop a5 pypi.org."),
    license = "BSD",
    keywords = "example documentation tutorial",
    url = "http://packages.python.org/an_example_pypi_project",
    packages=['an_example_pypi_project', 'tests'],
    long_description=read('README.md'),
    classifiers=[

    ],
)
