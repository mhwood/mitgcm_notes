# Testing MITgcm on a Computing Cluster

First, navigate to the directory where you'd like to perform your test and clone a fresh copy of MITgcm:
```
git clone https://github.com/MITgcm/MITgcm.git
```

In this example, we will use the `tutorial_barotropic_gyre` verification experiment to test the configuration on the cluster. Navigate to this directory and set up the code files to run with MPI:
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
module load intel mpich/intel hdf5/intel netcdf/intel
```
It is recommend that you purge all previously loaded modules before loading new ones:
```
module purge
```

### Choosing an opt file for the build
Next, design an optimization file for your cluster. There are a long list of options provided with MITgcm in the [MITgcm/tools/options](https://github.com/MITgcm/MITgcm/tree/master/tools/build_options) directory for a variety of systems.
For example, if you use the `intel` compiler above with `mpich`, your optfile may be 
```
linux_amd64_ifort
```
Note that you must also define an `MPI_INC_DIR` for your choice of MPI, e.g.:
```
export MPI_INC_DIR=/act/mpich/intel/include/
```

### Copmiling the model code
Next, move to the build directory and build your model with MPI using your optfile, e.g.:
```
../../../../tools/genmake2 -mpi -of ../../../../tools/build_options/linux_amd64_ifort -mo ../code
make depend
make
```

## Running the model
Now that the model has been built, test that it works using a job script for your cluster environment. For example, if your cluster uses slurm to manage jobs, your script called `job_tutorial_test` might looks like this:
```
#!/bin/bash
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=00:10:00
module purge
module load intel mpich/intel hdf5/intel netcdf/intel
ulimit -s unlimited
mpiexec -np 2 ./mitgcmuv
```

On slurm, submit this job to the queue as follows:
```
sbatch job_tutorial_test
```
