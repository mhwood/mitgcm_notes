# Compiling MITgcm on Spartan

At the time of this writing, MITgcm has been tested with the gnu compiler using the mpich MPI library along with the requisite netcdf and hdf5 libraries. Other software simiar to MITgcm has been compiled with other compilers (e.g. intel) and other MPI libraries (e.g. mvapich2, openmpi) so these *should* work similarly with MITgcm (if you'd like to help test, please get in touch!). To successfully compile MITgcm make sure the environment modules for each library (hdf5, netcdf, and one of the mpi libraries) are loaded before compiling.

## To compile MITgcm, use these steps: 
- First pick the compiler you want to use. I.E. Intel, PGI, or Gnu and version
- Load the appropriate module for compiler and libraries:

For gnu 6.3.0:
```
module load gnu/6.3.0 netcdf/gnu-6.3.0 mpich/gnu-6.3.0 hdf5/gnu-6.3.0
```

It's recommended that your purge previous modules before loading the above modules to avoid conflicts:
```
module purge all
```

