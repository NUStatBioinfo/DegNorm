from setuptools import setup

setup(
    name='DegNorm',
    version='0.0.1',
    packages=['degnorm'],
    entry_points={
        'console_scripts': ['degnorm=degnorm.__main__:main'],
    },
    package_data={
        'degnorm': ['resources/*'],
    },
    license='GNU',
    url='https://github.com/ffineis/DegNorm',
    author='Frank Fineis',
    author_email='frankfineis2022@u.northwestern.edu',
    description='DegNorm: the RNA-Seq read count normalization pipeline.'
)