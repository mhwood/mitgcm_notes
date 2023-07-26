# Compiling and Running MITgcm on Spartan

At the time of this writing, MITgcm has been tested with the gnu compiler using the mpich MPI library along with the requisite netcdf and hdf5 libraries. Other software simiar to MITgcm has been compiled with other compilers (e.g. intel) and other MPI libraries (e.g. mvapich2, openmpi) so these *should* work similarly with MITgcm (if you'd like to help test, please get in touch!). To successfully compile MITgcm make sure the environment modules for each library (hdf5, netcdf, and one of the mpi libraries) are loaded before compiling.

## To download MITgcm on Spartan, clone from the MITgcm github:
```
git clone https://github.com/MITgcm/MITgcm
```

## To compile MITgcm, use these steps: 
After cloning MITgcm, `cd` to the MITgcm directory and create a directory for your model. Within your models's directory, create a `code` and `build` directory. The `code` directory will hold modifications to the main MITgcm model code while the `build` directory will hold the compiled code and executable (mitgcmuv).

To begin, `cd` to the `build` directory.

It's recommended that your purge previous modules before loading the modules below to avoid conflicts:
```
module purge all
```

To compile:
- First pick the compiler you want to use. i.e. intel, PGI, or gnu and the version
- Load the appropriate module for compiler and libraries, for example:

For gnu 6.3.0:
```
module load gnu/6.3.0 netcdf/gnu-6.3.0 mpich/gnu-6.3.0 hdf5/gnu-6.3.0
```
Then, `export` a path to your MPI dir, e.g.:
```
export MPI_INC_DIR=/act/mpich/gnu-6.3.0/include/
```
Finally, use the following commands to compile, using the correct paths to `genmake2` and other MITgcm tools, e.g.:
```
../tools/genmake2 -mpi -of ../tools/build_options/linux_amd64_gfortran -mo ../code
make depend
make
```





