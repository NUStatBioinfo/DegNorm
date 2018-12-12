from setuptools import setup

setup(
    name='DegNorm',
    version='0.1.2',
    packages=['degnorm', 'degnorm.tests'],
    entry_points={
        'console_scripts': ['degnorm=degnorm.__main__:main',
                            'degnorm_mpi=degnorm.__main_mpi__:main',
                            'degnorm_test=degnorm.tests.__test__:main'],
    },
    package_data={
        'degnorm': ['resources/*', 'tests/data/*']
    },
    license='MIT',
    url='https://github.com/ffineis/DegNorm',
    author='Frank Fineis',
    author_email='frankfineis2022@u.northwestern.edu',
    description='DegNorm: the RNA-Seq read count normalization pipeline.',
    classifiers=[
        "Programming Language :: Python",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"]
)