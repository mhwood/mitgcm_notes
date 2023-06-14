# Testing MITgcm on a Computing Cluster

First, navigate to the directory where you'd like to perform your test and clone a fresh copy of MITgcm:
```
git clone https://github.com/MITgcm/MITgcm.git
```

In this example, we will use the `tutorial` verification experiment to test the configuration on the cluster. Navigate to this directory and set up the code files to run with MPI:
```
cd MITgcm/verification/tutorial_barotropic_gyre/
cd code
mv SIZE.h SIZE.h_no_mpi
mv SIZE.h_mpi SIZE.h
cd ..
```

## Building the model
Next, "build" (i.e. compile) the model. For this step, there are two pertinent items that need to be identified. 

### Choosing modules for the build
First, identify which compilers and associated modules are available on your cluster and load them accordingly. The command `module avail` will be helpful for this exercise.
For example, to use an intel compiler (ifort), your module loading might look like:
```
module load intel impi hdf5/intel netcdf/intel
```
If your compiler is `gnu`, your module load might look like:
```
module load gnu/9.3.0 mpich/gnu-9.3.0 hdf5/gnu-9.3.0 netcdf/gnu-9.3.0
```
It is recommend that you purge all previously loaded modules before loading new ones:
```
module purge
```

### Choosing an opt file for the build
Next, design an optimization file for your cluster. There are a long list of options provided with MITgcm in the [MITgcm/tools/options](https://github.com/MITgcm/MITgcm/tree/master/tools/build_options) directory for a variety of systems.
For example, if you use the `intel` compiler above with `impi`, your optfile may be 
```
linux_amd64_ifort+impi
```
If you use the `gnu` compiler with `mpich`, your optfile may be 
```
linux_amd64_ifort+gcc
```
Note that you must also define an `MPI_INC_DIR` for your choice of MPI, e.g.:
```
export MPI_INC_DIR=/act/mpich/intel/include/
```

## Building the model
Next, move to the build directory and build your model with MPI using your optfile, e.g.:
```
../../../../tools/genmake2 -mpi -of ../../../../tools/build_options/linux_amd64_ifort+gcc -mo ../code
make depend
make
```

