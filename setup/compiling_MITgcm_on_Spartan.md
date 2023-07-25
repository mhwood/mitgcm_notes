Compiling MITgcm on Spartan

All the prerequisites have already been compiled and installed on the cluster for each version of compiler. These prerequisites are jasper, hdf5, netcdf, and MPI libraries. Spartan has four different MPI libraries installed: impi (Intel MPI), mpich, mvapich2, and openmpi. They all should compile and work with WRF, but I have found that the mvapich2 library has the best performance. To successfully compile WRF you need to make sure that environment modules for each library (jasper, hdf5, netcdf, and one of the mpi libraries) are loaded before compiling WRF.
