from setuptools import setup

setup(
    name='DegNorm',
    version='0.0.1',
    packages=['degnorm'],
    entry_points = {
        'console_scripts': ['degnorm=degnorm.__main__:main'],
    },
    license='N/A'
)