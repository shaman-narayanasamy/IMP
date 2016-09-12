from setuptools import setup

setup(
    name='impy',
    version='0.1',
    py_modules=['impy'],
    install_requires=[
        'Click',
        'Path.py'
    ],
    entry_points='''
        [console_scripts]
        impy=impy:cli
    ''',
)
