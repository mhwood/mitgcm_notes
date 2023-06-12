# Installing MITgcm on MacOS

## Installing Open MPI

First download the latest version of OMPI [HERE](www.open-mpi.org/software/ompi/v4.1/). For MacOS, choose the `tar.gz` option.

Then, unzip the file and run these commandsL
```
commands
```

Next, add these lines to your `.bash_profile` file:
```
export MPI_INC_DIR="/Users/username/opt/usr/local/include"
export PATH="$PATH:/Users/username/opt/usr/local/bin"
export PKG_CONFIG_PATH="/Users/username/opt/usr/local/lib/pkgconfig"
export MPI_HOME="/Users/username/opt/usr/local/lib"
export TMPDIR="/tmp"
```
Replace `username` with your username on your machine.
