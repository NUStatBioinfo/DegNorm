# Installation

*This package is not yet on PyPi*. These instructions are for manual installation into a virtual environment.
DegNorm is only supported on *Nix platforms.

#### 1. Clone the DegNorm repository and `cd` into it.
```
git clone git@github.com:NUStatBioinfo/DegNorm.git
cd DegNorm
```

#### 2. Create a [conda](https://conda.io/docs/user-guide/tasks/manage-environments.html) virtual environment (or some other python3 virtual environment) and activate it:

    conda create -n degnorm python=3.6


#### 3. Run the `install` script

    ./install

# Testing

Check the successful installation of `degnorm` on your machine with the `degnorm_test` command. This runs all unit tests
and a minimal DegNorm pipeline run on a small batch of sample data.

By default, `degnorm_test` will clean up after itself by removing the temporary directory containing the output
of a full pipeline test run. If you would like to keep and inspect that directory, add the `--keep-output` flag:

    degnorm_test --keep-output
