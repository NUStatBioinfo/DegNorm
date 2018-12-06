# Installation

*This package is not yet on PyPi*. These instructions are for manual installation into a virtual environment.
DegNorm is only supported on *Nix platforms.

## Installing DegNorm in a conda environment

#### 1. Clone the DegNorm repository and `cd` into it.
```
git clone git@github.com:NUStatBioinfo/DegNorm.git
cd DegNorm
```

#### 2. Create a [conda](https://conda.io/docs/user-guide/tasks/manage-environments.html) virtual environment and activate it:

    conda create -n degnorm python=3.6


#### 3. Run the `install` script

    ./install

# Testing

Check the successful installation of `degnorm` on your machine with the `degnorm_test` command. This runs all unit tests
and a minimal DegNorm pipeline run on a small batch of sample data.

By default, `degnorm_test` will clean up after itself by removing the temporary directory containing the output
of a full pipeline test run. If you would like to keep and inspect that directory, add the `--keep-output` flag:

    degnorm_test --keep-output
    

## Requirements for `degnorm_mpi`

To make use of `degnorm_mpi`, you'll further need to install the [mpi4py](https://mpi4py.readthedocs.io/en/stable/index.html) package 
which requires the MPICH MPI library be installed and configured across your server system. Check that the `mpiexec` and/or 
`mpirun` commands are available to you at the command line. The `mpi4py` package will not be installed by default from the `.install` script.

If you're running DegNorm on a high performance computing environment (e.g. if you're at a university or research institute),
 it is most likely the case that MPI is already installed and available to you. If this is the case,
 then simply `pip install` the mpi4py package within the conda environment created in step 2:
    
    pip install mpi4py
 
Otherwise, you'll need to install [MPICH](https://en.wikipedia.org/wiki/MPICH).