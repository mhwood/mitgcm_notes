# Compiling MITgcm on Spartan

At the time of this writing, MITgcm has been tested with the gnu compiler using the mpich MPI library along with the requisite netcdf and hdf5 libraries. Other software simiar to MITgcm has been compiled with other compilers (e.g. intel) and other MPI libraries (e.g. mvapich2, openmpi) so these *should* work similarlu with MITgvm (if you'd like to help test, please get in touch!). To successfully compile MITgcm make sure the environment modules for each library (hdf5, netcdf, and one of the mpi libraries) are loaded before compiling.

## To compile MITgcm, use these steps: 
- First pick the compiler you want to use. I.E. Intel, PGI, or Gnu and version
- Load the appropriate module for compiler and libraries:
