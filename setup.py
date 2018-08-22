from setuptools import setup

setup(
    name='DegNorm',
    version='0.0.1',
    packages=['degnorm'],
    entry_points = {
        'console_scripts': ['degnorm=degnorm.__main__:main'],
    },
    package_data = {
        'degnorm': ['resources/*'],
    },
    license='N/A'
)